// lygemm -- optional runtime code generation.
//
// The compile time entry points only cover the shapes in LYGEMM_DISPATCH_LIST;
// everything else falls back to the blocked driver, which profiled 2.5x to
// 12.5x behind a kernel that is specialised for the exact shape (5x7x3: 1680
// cycles against 135).  This generates that kernel at run time: a straight
// line, fully unrolled AVX2 kernel with the displacements baked in, no loop,
// no branch, no mask built at run time, and the tiles picked from the shape.
//
//   C = A * B, row major and packed, alpha = 1 and beta = 0 only.
//
// Scope is deliberately the small end: up to 24 columns (six column registers)
// and a bounded amount of generated code.  Anything else returns null and the
// caller keeps its own path.  The kernels are cached per shape and are safe to
// call from several threads.
//
// SPDX-License-Identifier: MIT
#ifndef LYGEMM_JIT_HPP
#define LYGEMM_JIT_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include "lygemm/config.hpp"

#ifndef LYGEMM_USE_JIT
#  if defined(__linux__) && defined(__x86_64__) && LYGEMM_ARCH >= 2 && \
      LYGEMM_HAS_FMA
#    define LYGEMM_USE_JIT 1
#  else
#    define LYGEMM_USE_JIT 0
#  endif
#endif

#if LYGEMM_USE_JIT
#  include <sys/mman.h>
#  include <unistd.h>

#  ifndef LYGEMM_JIT_MAX_CODE
#    define LYGEMM_JIT_MAX_CODE 65536 /* bytes of generated code per kernel */
#  endif

#  include <atomic>
#  include <utility>

namespace lygemm {
namespace detail {
namespace jit {

/* ===================================================================== */
/* encoder: the six instruction forms a small GEMM kernel needs          */
/* ===================================================================== */

class emitter {
 public:
  explicit emitter(std::vector<std::uint8_t>& out) : out_(out) {}

  std::size_t here() const { return out_.size(); }
  void byte(std::uint8_t b) { out_.push_back(b); }
  void u32(std::uint32_t v) {
    for (int i = 0; i < 4; ++i) out_.push_back(std::uint8_t(v >> (8 * i)));
  }
  void patch_u32(std::size_t at, std::uint32_t v) {
    for (int i = 0; i < 4; ++i) out_[at + std::size_t(i)] = std::uint8_t(v >> (8 * i));
  }

 private:
  std::vector<std::uint8_t>& out_;
};

/* VEX3 prefix.  map: 1 = 0F, 2 = 0F38.  pp: 1 = the 66 prefix.  l: 1 for
 * 256-bit operands.  The R/X/B and vvvv fields are stored inverted. */
inline void vex3(emitter& e, int map, int pp, int w, int l, int vvvv, int reg,
                 int rm_reg, int rm_high) {
  const int r_bit = (reg >> 3) & 1;
  const int b_bit = rm_high;
  const int x_bit = 0; /* no SIB in the addressing modes used here */
  e.byte(0xC4);
  e.byte(std::uint8_t((((r_bit ^ 1) & 1) << 7) | (((x_bit ^ 1) & 1) << 6) |
                      (((b_bit ^ 1) & 1) << 5) | (map & 0x1F)));
  e.byte(std::uint8_t(((w & 1) << 7) | (((vvvv ^ 0xF) & 0xF) << 3) |
                      ((l & 1) << 2) | (pp & 3)));
  (void)rm_reg;
}

/* [base + disp]; the bases are rdi / rsi / rdx, so no SIB byte is needed.
 *
 * `tuple` is the disp8 scale of the instruction: 1 for VEX, where a disp8 is a
 * plain signed byte, and the tuple size for EVEX, where it is scaled (8 for a
 * broadcast, 64 for a full 512-bit operand).  Using the short form where it
 * fits cuts the code by about a fifth, which matters because these kernels are
 * straight line code long enough to be limited by instruction fetch. */
inline void modrm_mem_scaled(emitter& e, int reg, int base, std::int32_t disp,
                             int tuple) {
  if (tuple > 0 && disp >= 0 && disp % tuple == 0 && disp / tuple <= 127) {
    e.byte(std::uint8_t(0x40 | ((reg & 7) << 3) | (base & 7)));
    e.byte(std::uint8_t(disp / tuple));
    return;
  }
  e.byte(std::uint8_t(0x80 | ((reg & 7) << 3) | (base & 7)));
  e.u32(static_cast<std::uint32_t>(disp));
}

inline void modrm_mem(emitter& e, int reg, int base, std::int32_t disp) {
  modrm_mem_scaled(e, reg, base, disp, 1);
}

inline void modrm_reg(emitter& e, int reg, int rm) {
  e.byte(std::uint8_t(0xC0 | ((reg & 7) << 3) | (rm & 7)));
}

/* three register operands */
inline void vreg3(emitter& e, int opcode, int map, int pp, int w, int dst,
                  int src1, int src2) {
  vex3(e, map, pp, w, 1, src1, dst, src2, (src2 >> 3) & 1);
  e.byte(std::uint8_t(opcode));
  modrm_reg(e, dst, src2);
}

inline void zero(emitter& e, int r) { vreg3(e, 0x57, 1, 1, 0, r, r, r); }

/* dst = dst * a + b  (AVX2 FMA cannot take a memory operand unaligned, so the
 * operands are always registers) */
inline void fma231(emitter& e, int dst, int a, int b) {
  vreg3(e, 0xB8, 2, 1, 1, dst, a, b);
}

inline void load(emitter& e, int dst, int base, std::int32_t disp) {
  vex3(e, 1, 1, 0, 1, 0, dst, 0, 0);
  e.byte(0x10);
  modrm_mem(e, dst, base, disp);
}

inline void store(emitter& e, int base, std::int32_t disp, int src) {
  vex3(e, 1, 1, 0, 1, 0, src, 0, 0);
  e.byte(0x11);
  modrm_mem(e, src, base, disp);
}

inline void broadcast(emitter& e, int dst, int base, std::int32_t disp) {
  vex3(e, 2, 1, 0, 1, 0, dst, 0, 0);
  e.byte(0x19);
  modrm_mem(e, dst, base, disp);
}

inline void mask_load(emitter& e, int dst, int mask, int base,
                      std::int32_t disp) {
  vex3(e, 2, 1, 0, 1, mask, dst, 0, 0);
  e.byte(0x2D);
  modrm_mem(e, dst, base, disp);
}

/* VPMASKMOVPD m256, ymm1, ymm2: the *data* is in ModRM.reg and the mask in
 * VEX.vvvv, while the load form is the other way round.  Getting this wrong
 * makes the stored lanes depend on the sign bits of the accumulator, which
 * writes almost nothing and only for negative results. */
inline void mask_store(emitter& e, int base, std::int32_t disp, int mask,
                       int src) {
  vex3(e, 2, 1, 0, 1, mask, src, 0, 0);
  e.byte(0x2F);
  modrm_mem(e, src, base, disp);
}

/* rip-relative load, used once for the tail mask */
inline void load_rip(emitter& e, int dst, std::size_t& patch_at) {
  vex3(e, 1, 1, 0, 1, 0, dst, 0, 0);
  e.byte(0x10);
  e.byte(std::uint8_t(0x05 | ((dst & 7) << 3)));
  patch_at = e.here();
  e.u32(0);
}

inline void vzeroupper_(emitter& e) {
  e.byte(0xC5);
  e.byte(0xF8);
  e.byte(0x77);
}

inline void ret_(emitter& e) { e.byte(0xC3); }

#if LYGEMM_ARCH >= 3
/* ---------------------------------------------------------------- AVX-512
 *
 * Four byte EVEX prefix.  Vector length is fixed at 512 bits, the mask comes
 * from a k register, and `bcst` asks for an embedded broadcast, which lets the
 * A operand of an FMA come straight from memory -- one instruction per (row,
 * k) instead of a load, a broadcast and an FMA.
 *
 * Registers 16..31 are reached through R' (bit 4 of reg), X (bit 4 of a
 * register rm) and V' (bit 4 of vvvv).  Leaving V' set costs nothing below 16
 * registers and silently makes every wider register alias to reg & 15 above it,
 * which is how the 512-bit generator first produced wrong results. */
inline void evex(emitter& e, int mm, int pp, int w, int vvvv, int bcst, int z,
                 int mask, int reg, int rm_high, int rm_high4) {
  const int r_bit = (reg >> 3) & 1;
  const int r4_bit = (reg >> 4) & 1;
  const int v4_bit = (vvvv >> 4) & 1;
  e.byte(0x62);
  e.byte(std::uint8_t((((r_bit ^ 1) & 1) << 7) | (((rm_high4 ^ 1) & 1) << 6) |
                      (((rm_high ^ 1) & 1) << 5) | (((r4_bit ^ 1) & 1) << 4) |
                      (mm & 3)));
  e.byte(std::uint8_t(((w & 1) << 7) | (((vvvv ^ 0xF) & 0xF) << 3) | (1 << 2) |
                      (pp & 3)));
  e.byte(std::uint8_t(((z & 1) << 7) | (1 << 6) /* L'L = 512 bits */
                      | ((bcst & 1) << 4) | (((v4_bit ^ 1) & 1) << 3) |
                      (mask & 7)));
}

inline void zmm_reg3(emitter& e, int opcode, int map, int pp, int w, int dst,
                     int src1, int src2) {
  evex(e, map, pp, w, src1, 0, 0, 0, dst, (src2 >> 3) & 1, (src2 >> 4) & 1);
  e.byte(std::uint8_t(opcode));
  modrm_reg(e, dst, src2);
}

inline void zmm_zero(emitter& e, int r) { zmm_reg3(e, 0x57, 1, 0, 0, r, r, r); }

/* dst += src1 * [base + disp], the A operand broadcast from memory.  Reading
 * the same broadcast once per column register is only worth it for a single
 * vector; from two up the broadcast is hoisted into a register instead. */
inline void zmm_fma_bcst(emitter& e, int dst, int src1, int base,
                         std::int32_t disp, int mask) {
  evex(e, 2, 1, 1, src1, 1 /* broadcast */, 0, mask, dst, 0, 0);
  e.byte(0xB8);
  modrm_mem_scaled(e, dst, base, disp, 8);
}

inline void zmm_fma_reg(emitter& e, int dst, int src1, int src2, int mask) {
  evex(e, 2, 1, 1, src1, 0, 0, mask, dst, (src2 >> 3) & 1, (src2 >> 4) & 1);
  e.byte(0xB8);
  modrm_reg(e, dst, src2);
}

/* dst = [base + disp]; the 512-bit forms of the moves carry W = 1 (the FMA
 * uses W = 1 too, but for a different reason, and vxorps uses W = 0). */
inline void zmm_bcast(emitter& e, int dst, int base, std::int32_t disp) {
  evex(e, 2, 1, 1, 0, 0, 0, 0, dst, 0, 0);
  e.byte(0x19);
  modrm_mem_scaled(e, dst, base, disp, 8);
}

inline void zmm_load(emitter& e, int dst, int base, std::int32_t disp,
                     int mask) {
  evex(e, 1, 1, 1, 0, 0, mask != 0 ? 1 : 0, mask, dst, 0, 0);
  e.byte(0x10);
  modrm_mem_scaled(e, dst, base, disp, 64);
}

inline void zmm_store(emitter& e, int base, std::int32_t disp, int src,
                      int mask) {
  evex(e, 1, 1, 1, 0, 0, 0, mask, src, 0, 0);
  e.byte(0x11);
  modrm_mem_scaled(e, src, base, disp, 64);
}

/* k1 = [rip + disp32], a 16 bit mask constant emitted after the code */
inline void kmovw_rip(emitter& e, int kreg, std::size_t& patch_at) {
  e.byte(0xC5);
  e.byte(0xF8); /* VEX.128.0F.W0, no vvvv */
  e.byte(0x90); /* kmovw */
  e.byte(std::uint8_t(0x05 | ((kreg & 7) << 3)));
  patch_at = e.here();
  e.u32(0);
}
#endif /* LYGEMM_ARCH >= 3 */

/* ===================================================================== */
/* kernel generation                                                     */
/* ===================================================================== */

/* Registers: rdi = A, rsi = B, rdx = C (SysV).  ymm0..15 hold the whole tile,
 * so the kernel is a leaf with no prologue and no stack. */
enum { reg_a = 7, reg_b = 6, reg_c = 2, max_ymm = 16 };

inline std::int32_t disp(std::size_t elem) {
  return static_cast<std::int32_t>(elem * sizeof(double));
}

#if LYGEMM_ARCH >= 3
/* 512-bit kernel.  The A operand is broadcast straight out of memory inside
 * the FMA, so a (row, k) step costs one instruction rather than the 256-bit
 * path's load plus broadcast plus FMA, and the column tail rides in a k mask.
 * Thirty two registers also buy a much taller tile: six column registers still
 * leave four rows, where the 256-bit path fits one. */
inline std::size_t generate_avx512(std::vector<std::uint8_t>& out,
                                   std::size_t m, std::size_t n,
                                   std::size_t k) {
  const std::size_t nv = (n + 7) / 8;      /* column registers */
  const std::size_t tail = n - (nv - 1) * 8; /* live lanes in the last */
  if (nv > 8) return 0;
  const bool need_mask = (tail != 8);
  const std::size_t bcast_reg = (nv > 1) ? 1u : 0u; /* A broadcast in a register */
  std::size_t rows = (32 - nv - bcast_reg) / nv;
  if (rows > 16) rows = 16;
  if (rows > m) rows = m;
  if (rows == 0) return 0;
  const std::size_t tiles = (m + rows - 1) / rows;
  if (tiles * k * (rows + 1) * nv * 8 > LYGEMM_JIT_MAX_CODE) return 0;

  emitter e(out);
  std::size_t mask_patch = 0;
  if (need_mask) kmovw_rip(e, 1, mask_patch);

  for (std::size_t i0 = 0; i0 < m; i0 += rows) {
    const std::size_t r = (i0 + rows <= m) ? rows : (m - i0);
    for (std::size_t i = 0; i < r; ++i)
      for (std::size_t v = 0; v < nv; ++v) zmm_zero(e, int(i * nv + v));

    for (std::size_t kk = 0; kk < k; ++kk) {
      for (std::size_t v = 0; v < nv; ++v) {
        const int dst = int(r * nv + v);
        const int msk = (need_mask && v + 1 == nv) ? 1 : 0;
        zmm_load(e, dst, reg_b, disp(kk * n + v * 8), msk);
      }
      for (std::size_t i = 0; i < r; ++i) {
        const std::int32_t a = disp((i0 + i) * k + kk);
        if (nv == 1) {
          zmm_fma_bcst(e, int(i * nv), int(r * nv), reg_a, a, need_mask ? 1 : 0);
        } else {
          const int bc = int(r * nv + nv);
          zmm_bcast(e, bc, reg_a, a);
          for (std::size_t v = 0; v < nv; ++v)
            zmm_fma_reg(e, int(i * nv + v), bc, int(r * nv + v),
                        (need_mask && v + 1 == nv) ? 1 : 0);
        }
      }
    }

    for (std::size_t i = 0; i < r; ++i)
      for (std::size_t v = 0; v < nv; ++v) {
        const int msk = (need_mask && v + 1 == nv) ? 1 : 0;
        zmm_store(e, reg_c, disp((i0 + i) * n + v * 8), int(i * nv + v), msk);
      }
  }

  vzeroupper_(e);
  ret_(e);

  if (need_mask) {
    while ((out.size() & 15U) != 0) e.byte(0x90);
    const std::size_t at = out.size();
    const std::uint16_t bits = std::uint16_t((1u << tail) - 1u);
    e.byte(std::uint8_t(bits & 0xFF));
    e.byte(std::uint8_t(bits >> 8));
    const std::int64_t rel = std::int64_t(at) - std::int64_t(mask_patch + 4);
    e.patch_u32(mask_patch, std::uint32_t(rel));
  }
  return out.size();
}
#endif /* LYGEMM_ARCH >= 3 */

/* Emits C = A*B for a packed row major m x n x k.  Returns the code size, or 0
 * when the shape is outside the supported range. */
inline std::size_t generate(std::vector<std::uint8_t>& out, std::size_t m,
                            std::size_t n, std::size_t k) {
  out.clear();
  if (m == 0 || n == 0 || k == 0) return 0;
  if (m > 0xFFFF || n > 0xFFFF || k > 0xFFFF) return 0;
#if LYGEMM_ARCH >= 3
  /* five columns and up fill a 512-bit register well enough to beat two
   * 256-bit ones; narrower shapes stay on the 256-bit path */
  if (n >= 5) return generate_avx512(out, m, n, k);
#endif

  const std::size_t nv = (n + 3) / 4;      /* column registers */
  const std::size_t tail = n - (nv - 1) * 4; /* live lanes in the last one */
  if (nv > 6) return 0;
  const bool need_mask = (tail != 4);

  /* accumulators + B row + broadcast + mask */
  const std::size_t fixed = nv + 1 + (need_mask ? 1 : 0);
  if (fixed + nv > max_ymm) return 0;
  std::size_t rows = (max_ymm - fixed) / nv;
  if (rows > m) rows = m;
  if (rows == 0) return 0;

  const std::size_t tiles = (m + rows - 1) / rows;
  /* a rough bound on the code before generating it */
  const std::size_t per_k = nv + rows + rows * nv;
  if (tiles * k * per_k * 8 > LYGEMM_JIT_MAX_CODE) return 0;

  emitter e(out);
  const int mask_reg = static_cast<int>(rows * nv + nv + 1);
  std::size_t mask_patch = 0;
  if (need_mask) load_rip(e, mask_reg, mask_patch);

  for (std::size_t i0 = 0; i0 < m; i0 += rows) {
    const std::size_t r = (i0 + rows <= m) ? rows : (m - i0);

    for (std::size_t i = 0; i < r; ++i)
      for (std::size_t v = 0; v < nv; ++v)
        zero(e, static_cast<int>(i * nv + v));

    for (std::size_t kk = 0; kk < k; ++kk) {
      /* the B row stays in registers for the whole tile */
      for (std::size_t v = 0; v < nv; ++v) {
        const int dst = static_cast<int>(r * nv + v);
        const std::int32_t d = disp(kk * n + v * 4);
        if (need_mask && v + 1 == nv)
          mask_load(e, dst, mask_reg, reg_b, d);
        else
          load(e, dst, reg_b, d);
      }
      for (std::size_t i = 0; i < r; ++i) {
        const int bc = static_cast<int>(r * nv + nv);
        broadcast(e, bc, reg_a, disp((i0 + i) * k + kk));
        for (std::size_t v = 0; v < nv; ++v) {
          const int acc = static_cast<int>(i * nv + v);
          fma231(e, acc, bc, static_cast<int>(r * nv + v));
        }
      }
    }

    for (std::size_t i = 0; i < r; ++i)
      for (std::size_t v = 0; v < nv; ++v) {
        const int acc = static_cast<int>(i * nv + v);
        const std::int32_t d = disp((i0 + i) * n + v * 4);
        if (need_mask && v + 1 == nv)
          mask_store(e, reg_c, d, mask_reg, acc);
        else
          store(e, reg_c, d, acc);
      }
  }

  vzeroupper_(e);
  ret_(e);

  if (need_mask) {
    /* the constant follows the code, 32 byte aligned */
    while ((out.size() & 31U) != 0) e.byte(0x90);
    const std::size_t at = out.size();
    for (std::size_t lane = 0; lane < 4; ++lane) {
      const std::uint64_t v = (lane < tail) ? ~std::uint64_t(0) : 0;
      for (int i = 0; i < 8; ++i)
        e.byte(std::uint8_t(v >> (8 * i)));
    }
    const std::int64_t rel =
        static_cast<std::int64_t>(at) - static_cast<std::int64_t>(mask_patch + 4);
    e.patch_u32(mask_patch, static_cast<std::uint32_t>(rel));
  }
  return out.size();
}

/* ===================================================================== */
/* cache                                                                 */
/* ===================================================================== */

using kernel_fn = void (*)(const double*, const double*, double*);

inline std::uint64_t pack_key(std::size_t m, std::size_t n, std::size_t k) {
  return (std::uint64_t(m) << 40) | (std::uint64_t(n) << 20) | std::uint64_t(k);
}

struct slot {
  std::atomic<std::uint64_t> key{0};
  std::atomic<kernel_fn> fn{nullptr};
};

struct cache {
  /* Open addressing, never evicted: a program that asks for more distinct
   * shapes than this falls back to the C++ path rather than thrashing.  A test
   * suite walking every shape from 1x1x1 to 17x17x17 fits. */
  static constexpr std::size_t size = 8192;
  slot slots[size];
  std::atomic<bool> full{false};
};

/* Namespace scope so that the lookup does not pay for a function local static's
 * initialisation guard; the atomics make it a constant initialiser. */
inline cache the_cache;

/* Copies the code into its own executable mapping. */
inline kernel_fn materialise(const std::vector<std::uint8_t>& code) {
  if (code.empty()) return nullptr;
  const std::size_t page = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
  const std::size_t bytes = ((code.size() + page - 1) / page) * page;
  void* mem = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (mem == MAP_FAILED) return nullptr;
  std::memcpy(mem, code.data(), code.size());
  __builtin___clear_cache(static_cast<char*>(mem),
                          static_cast<char*>(mem) + code.size());
  if (mprotect(mem, bytes, PROT_READ | PROT_EXEC) != 0) {
    munmap(mem, bytes);
    return nullptr;
  }
  return reinterpret_cast<kernel_fn>(mem);
}

inline std::size_t hash_key(std::uint64_t key) {
  return static_cast<std::size_t>(
      (key * 0x9E3779B97F4A7C15ull) >> (64 - 9));
}

/* Generates on first use and caches.  Safe from several threads: two threads
 * may both build the same kernel, one wins the slot and the other's code is
 * dropped. */
inline kernel_fn kernel(std::size_t m, std::size_t n, std::size_t k) {
  const std::uint64_t key = pack_key(m, n, k);
  cache& c = the_cache;
  std::size_t i = hash_key(key) & (cache::size - 1);
  for (std::size_t probe = 0; probe < cache::size; ++probe) {
    const std::uint64_t at = c.slots[i].key.load(std::memory_order_acquire);
    if (at == key) return c.slots[i].fn.load(std::memory_order_acquire);
    if (at == 0) {
      kernel_fn mine = nullptr;
      {
        static thread_local std::vector<std::uint8_t> code;
        if (generate(code, m, n, k) == 0) return nullptr;
        mine = materialise(code);
        if (mine == nullptr) return nullptr;
      }
      std::uint64_t expected = 0;
      if (c.slots[i].key.compare_exchange_strong(expected, key,
                                                 std::memory_order_acq_rel)) {
        c.slots[i].fn.store(mine, std::memory_order_release);
        return mine;
      }
      
      if (expected == key) {
        /* another thread published first: use theirs */
        kernel_fn other = c.slots[i].fn.load(std::memory_order_acquire);
        if (other != nullptr) return other;
      }
      /* the slot was taken by a different shape: probe on */
    }
    i = (i + 1) & (cache::size - 1);
  }
  return nullptr;
}

inline bool enabled() { return true; }

/* Typed entry point; the generated kernels are for double only for now. */
template <class T>
struct entry {
  using fn = void (*)(const T*, const T*, T*);
  static fn get(std::size_t, std::size_t, std::size_t) { return nullptr; }
};

template <>
struct entry<double> {
  using fn = void (*)(const double*, const double*, double*);
  static fn get(std::size_t m, std::size_t n, std::size_t k) {
    return kernel(m, n, k);
  }
};

} // namespace jit
} // namespace detail
} // namespace lygemm

#else /* LYGEMM_USE_JIT */

namespace lygemm {
namespace detail {
namespace jit {

inline bool enabled() { return false; }

template <class T>
struct entry {
  using fn = void (*)(const T*, const T*, T*);
  static fn get(std::size_t, std::size_t, std::size_t) { return nullptr; }
};

} // namespace jit
} // namespace detail
} // namespace lygemm

#endif /* LYGEMM_USE_JIT */

#endif /* LYGEMM_JIT_HPP */
