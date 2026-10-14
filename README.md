# lygemm

Header-only small matrix multiplication for x86-64, C++17. No build step, no
dependencies: `-Iinclude` and multiply.

The kernel is organised the way tile DSLs are, with the differences a CPU
makes: a GEMM is a grid of output tiles, every tile owns a `fragment` of
registers, and the K loop feeds panels of A and B into it.

```cpp
#include <lygemm.hpp>

double A[16] = {...}, B[16] = {...}, C[16];
lygemm::gemm<4, 4, 4>(A, B, C);              // compile time shape
lygemm::gemm(m, n, k, A, B, C);              // runtime shape
lygemm::gemm_batch<8, 8, 8>(count, A, B, C); // `count` independent products
```

## What it is for

Small matrices, the sizes that spend more time in call overhead than in
arithmetic. `${}^3\sqrt{MNK} \le 64$` is the usual boundary definition
([LIBXSMM](https://libxsmm.github.io/libxsmm/documentation/)).

| Size | Where it comes from |
| --- | --- |
| 1x1 … 8x8 | fixed-size Kalman filtering, streaming sensor fusion ([fast-kalman](https://fast-kalman.dev/)) |
| 4x4 | rigid body transforms (SE(3)), graphics, EKF state covariance ([lecture notes](https://cecas.clemson.edu/~ahoover/ece854/lecture-notes/lecture-ekf.pdf)) |
| 3x3, 5x5, 7x7 | convolution kernels; `K = C·k·k` after im2col ([GEMM-based convolution](https://arxiv.org/pdf/1709.03395)) |
| 4-8 x 8-16 | deep learning micro kernels ([XNNPACK](https://chromium.googlesource.com/external/github.com/google/XNNPACK/) 8x8, 1x16) |
| 6x6, 12x12, 24x24 | 6-DOF adjoints, inertial state, covariance blocks |

Above roughly 16-24 the general libraries win and this one does not try: hand
the work to BLAS (see `LYGEMM_USE_BLAS` below).

## API

All matrices are row major. `C` must not overlap `A` or `B`. Shapes must be
non zero. `T` is `double` or `float`.

```cpp
namespace lygemm {

// compile time shapes, fully specialised
template <std::size_t M, std::size_t N, std::size_t K, class T = double>
void gemm(const T* A, const T* B, T* C);
template <std::size_t M, std::size_t N, std::size_t K, class T = double>
void gemm(T alpha, const T* A, const T* B, T beta, T* C);

// runtime shapes
template <class T>
void gemm(std::size_t m, std::size_t n, std::size_t k, const T* A, const T* B,
          T* C, T alpha = 1, T beta = 0);

// batches, packed or strided; a stride of 0 broadcasts that operand
template <std::size_t M, std::size_t N, std::size_t K, class T = double>
void gemm_batch(std::size_t batch, const T* A, const T* B, T* C, T alpha = 1, T beta = 0);
template <std::size_t M, std::size_t N, std::size_t K, class T = double>
void gemm_batch(std::size_t batch, const T* A, std::ptrdiff_t sa, const T* B,
                std::ptrdiff_t sb, T* C, std::ptrdiff_t sc, T alpha = 1, T beta = 0);
template <class T>
void gemm_batch(std::size_t batch, std::size_t m, std::size_t n, std::size_t k, ...);

// one operand shared by the whole batch: C[b] = alpha * A * B[b] + beta * C[b]
template <std::size_t M, std::size_t N, std::size_t K, class T = double>
void gemm_batch_shared_a(std::size_t batch, const T* A, const T* B, std::ptrdiff_t sb,
                         T* C, std::ptrdiff_t sc, T alpha = 1, T beta = 0);
template <std::size_t M, std::size_t N, std::size_t K, class T = double>
void gemm_batch_shared_b(std::size_t batch, const T* A, std::ptrdiff_t sa, const T* B,
                         T* C, std::ptrdiff_t sc, T alpha = 1, T beta = 0);

} // namespace lygemm
```

`beta == 0` never reads `C`, so an uninitialised output is fine.

### Getting a kernel to call in a loop

For a shape used repeatedly, fetch the generated kernel once and call it
directly.  That skips the per call dispatch and lookup, and is the fastest way
to use this library:

```cpp
using fn = lygemm::kernel_fn<double>;
fn k4 = lygemm::kernel<double>(4, 4, 4);   // null if unavailable, see below
if (k4) for (int i = 0; i < n; ++i) k4(A + i * 16, B + i * 16, C + i * 16);
else    for (int i = 0; i < n; ++i) lygemm::gemm<4, 4, 4>(A + i * 16, B + i * 16, C + i * 16);
```

The generated kernels cover `double`, `alpha == 1 && beta == 0`, packed row
major operands, up to eight column registers and a bounded amount of code;
`lygemm::kernel` returns null for anything else, and the compile time entry
points and `lygemm::gemm` handle the rest.  A kernel is valid for the life of
the process.

Sharing one operand across the batch hoists the tile grid out of the batch
loop, so it is worth expressing. Do not expect it to be faster on these sizes
though: a small matrix is L1 resident, so the panel loads it saves were already
cheap. Measured, `gemm_batch_shared_a` and the plain strided form with a zero
stride come out within noise of each other (8x8x8: 98.6 ns against 97.8 ns per
item).

### Tile level API

The pieces the kernels are built from are public, for when the default grid is
not what you want:

```cpp
lygemm::tile::fragment<double, 4, 8> acc;   // register tile of the output
acc.clear();
lygemm::tile::mma<8>(A_panel, lda, B_panel, ldb, acc);  // acc += A_panel * B_panel
acc.store(C_tile, ldc, alpha, beta);
```

## How it is built

The compile time path is fully specialised:

```
gemm<M,N,K>(A,B,C)
  grid over output tiles (runtime loop, tile shape is the template argument)
    fragment<BM,BN> acc; acc.clear();
    for each K panel of BK: acc += A_panel * B_panel
    acc.store(C_tile, ldc, alpha, beta)      // C = alpha*acc + beta*C
```

1. **Register blocking.** The tile is three column registers wide (`3*lanes`
   columns, 12 for doubles) and as many rows as 12 accumulators allow. Covering
   all of M in one tile where the registers allow is worth 33% on 6x6x6, 36% on
   5x5x5 and 32% on 8x1x8 against stopping at four rows, because every extra
   row tile re-reads the B panel and repeats the fragment epilogue.
2. **Runtime grid, compile time tiles.** Unrolling the grid instead multiplies
   the code by the number of tiles; measured on this machine, 16x16x16 is 4x
   faster with a runtime grid than with the unrolled one.
3. **Masked column tail.** When N is not a multiple of the vector width the
   last column register is loaded and stored with a mask. A scalar tail looks
   cheaper, but the vectoriser rewrites it into shuffles: 3x3x3 went from 53
   instructions (24 shuffles, 9 FMAs) to 12.
4. **Always inlined unrolling.** `static_for` and the tile bodies are marked
   `always_inline`. If the compiler is allowed to outline one of them the
   fragment is passed by reference and its register tile is demoted to the
   stack, which costs several times the kernel.
5. **Runtime shapes** first look for a shape in `LYGEMM_DISPATCH_LIST` (squares
   1..8, 12, 16) and otherwise run a blocked kernel that stages partial panels
   through zero padded scratch, so no access is ever masked or guarded.
6. **Per shape vector width.** A backend wider than the column count is a loss:
   four columns leave half of an AVX-512 register idle in every FMA, and the
   parts that have AVX-512 downclock for 512-bit work. `shape_ops<T, N>` picks
   the widest backend the columns fill, so on a Skylake-SP a 4x4x4 goes through
   the 256-bit path in both builds: 20.5 ns, against 43.9 ns before.
7. **A kernel generated for the exact shape (optional).** `LYGEMM_USE_JIT` is
   on by default on x86-64 with AVX2 and FMA.  The first call for a shape emits
   a straight line kernel with every displacement baked in -- no loop, no
   branch, no mask built at run time -- caches it, and calls it through a
   pointer.  Code is generated at the build's own width: 512-bit EVEX kernels
   when compiled for AVX-512 (where the A operand is broadcast straight out of
   memory by the FMA itself, one instruction per (row, k)), 256-bit VEX below
   that.  On a Xeon Gold 6132 that is 89 instructions for an 8x8x8 where the
   compile time path needs 442 and libxsmm's generated kernel 143.

   `lygemm::kernel<T>(m, n, k)` hands out the kernel so a fixed shape can be
   called without any dispatch at all; that is the fastest way to use the
   library, and the way the benchmarks compare it against libxsmm.

   The limits are deliberate: `double`, `alpha == 1 && beta == 0`, packed row
   major operands, up to eight column registers, 256-bit and 512-bit vectors,
   and a bounded amount of code per kernel.  Everything else falls through to
   the paths above.  Kernels are cached per shape in a fixed table; a program
   that asks for more distinct shapes than it holds keeps working on the C++
   paths.

   Getting an encoder right is not obvious, and the two mistakes made here are
   worth recording: the store form of `vmaskmovpd` keeps the data in ModRM.reg
   and the mask in VEX.vvvv while the load form is the other way round, and in
   EVEX the bits that extend a register past 15 (R', X, V') are not optional --
   leaving V' set makes every register above 15 alias to `reg & 15`.  Both
   produced wrong results rather than faults.  The test suite therefore calls
   `detail::jit::entry<T>::get` directly on every shape from 1x1x1 to 10x10x10
   and compares against the reference on guard page buffers, so a silently
   declined shape cannot hide a broken encoder.

### Why intrinsics and not inline assembly

The two places the compiler got in the way are fixed structurally, in points 2
and 4 above, and what is left is close to the metal. Against the FMA throughput
bound of 0.5 cycles per 4-wide FMA, the kernels run at

| shape | 2x2x2 | 4x4x4 | 8x8x8 | 4x4x16 | 12x12x12 | 16x16x16 | 24x24x24 | 64x64x64 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cycles per FMA | 0.6 | 0.5 | 0.5 | 0.5-0.7 | 0.5-0.7 | 0.6 | 0.6 | 0.6 |

with no calls and no spills in the generated code, so there is no free slot
left for hand scheduling to win, and intrinsics map one to one onto
instructions anyway. Assembly would cost the portability across GCC, clang,
the three vector backends and the scalar fallback that the header-only form
exists for. It stays an option for a kernel that needs an instruction order the
register allocator will not produce; the two cases where the compiler actually
did fail are recorded above, so the next person can judge.

Instruction set: AVX-512 (8 doubles / 16 floats) when compiled for it, else
AVX2 (4/8), else portable scalar. Selected at compile time from the flags, so
`-mavx2 -mfma` or `-march=native` is all that is needed.

## Performance

Two machines, `-O3`, single core. The ceiling is
[cpufp](https://github.com/pigirons/cpufp)'s f64 FMA figure.

| | laptop (Arrow Lake-H) | server (Xeon Gold 6132) |
| --- | --- | --- |
| vector width | AVX2, 4 doubles | AVX-512, 8 doubles |
| f64 FMA ceiling, one core | 58.7 GFLOP/s | 53.8 GFLOP/s |

### Against libxsmm

Cycles per call on the server, one core, AVX-512 build, all three entry points
of this library against libxsmm's kernel for the same shape. Nothing else is
close to libxsmm on this range, so it is the one to beat.

| shape | compile time `gemm<M,N,K>` | `gemm(m,n,k,...)` | `kernel<double>(m,n,k)` | libxsmm |
| --- | --- | --- | --- | --- |
| 4x4x4 | 34 | 63 | **47** | 59 |
| 6x6x6 | 176 | 75 | **48** | 59 |
| 8x8x8 | 240 | 86 | **59** | 76 |
| 12x12x12 | 821 | 216 | **192** | 202 |
| 16x16x16 | 1482 | 413 | 392 | 352 |
| 24x24x24 | 4486 | 1198 | 1171 | 1049 |

Through `lygemm::kernel` this library is **level with libxsmm or ahead of it** from
4x4x4 to 12x12x12 (0.77x to 0.95x of libxsmm's cycles) and within 11% at
16x16x16 and 24x24x24.  Through `lygemm::gemm(m, n, k, ...)` it pays about 25
to 60 cycles for the dispatch and lookup on every call, and the compile time
entry points are for the shapes you know at build time.

Why the generated kernel wins where it does: the FMA takes its A operand as an
embedded broadcast straight out of memory, so a single column register costs one
instruction per (row, k) where the compile time path needs a load, a broadcast
and an FMA.  From two column registers up the broadcast is hoisted into a
register instead, because reading the same element three times costs more than
the extra instruction (24x24x24: 1171 cycles against 2232 with the broadcast
left embedded, and 797 L1 loads against 1949).  In instruction count the 8x8x8
is 98 against libxsmm's 143.

### Against the libraries, on the laptop

* **Against a naive triple loop**, same run, on the server: 2x2x2 2.9x, 3x3x3
  2.9x, 4x4x4 2.4x, 4x4x16 2.5x faster.
* **CBLAS**: a 4x4x4 costs 39 ns through the CBLAS call, against 2.0 ns for
  this library, and that call overhead is what the small shapes actually pay.
* From 12x12x12 the general libraries win. That crossover is the point of this
  library.
* The `naive` rows of the laptop table are auto-vectorised, because their sizes
  are compile time constants, so read the server rows for that comparison.

### Reading these numbers

* **Do not quote the nanoseconds without re-measuring.** The laptop is a VM
  whose clock moves by about 3x during a long run; the server idles at 1.2 GHz
  and turboes only under load. Absolute times from either moved by 2-3x between
  runs while the ratios above held, and even `rdtsc` based cycle counts move
  with the core, since the TSC ticks at a fixed reference rate.
  `bench/bench_gemm.cpp` rotates eight input buffer sets so no product can be
  hoisted out of the timing loop; the perf counters above come from `perf stat`
  on a pinned core.
* **Shapes narrower than the vector width do not benefit.** A 2x2x2 uses one
  quarter of each 256-bit register, so a scalar loop keeps up with it.
* **AVX2 can beat AVX-512 on the small end**, which is why the backend is
  picked per shape (point 6) rather than by the build flags.

## Building

There is nothing to build. To use the library, add `include/` to the include
path.

```sh
make test        # correctness suite
make test-asan   # the same, under ASan + UBSan and guard pages
make bench       # comparison benchmark, needs `bench/fetch_deps.sh` once
```

With CMake, either `add_subdirectory(lygemm)` and link `lygemm::lygemm`, or

```sh
cmake -S . -B build && cmake --build build && ctest --test-dir build
```

Verified on an Arrow Lake-H laptop (AVX2, GCC 13, clang 18, C++17 and C++20),
on a Xeon Gold 6132 (AVX-512, GCC 11), and with the scalar backend
(`-DLYGEMM_ARCH=0`), float and double. `LYGEMM_USE_BLAS` is built and run
against OpenBLAS.

The correctness suite covers every shape from 1x1x1 to 17x17x17, the larger
shapes up to 64x64x64, alpha/beta combinations, float and double, batched and
shared-operand batching, and the tile API. Buffers are `mmap`ed so that their
last element ends exactly at a `PROT_NONE` page: an out of bounds access, by a
SIMD load or otherwise, faults immediately instead of silently corrupting
memory.

## Configuration

| macro | default | meaning |
| --- | --- | --- |
| `LYGEMM_ARCH` | from the flags | 0 scalar, 2 AVX2, 3 AVX-512 |
| `LYGEMM_MAX_UNROLL_K` | 32 | K panels unrolled up to this, looped above |
| `LYGEMM_DISPATCH_LIST` | squares 1..8, 12, 16 | shapes the runtime entry point specialises |
| `LYGEMM_USE_BLAS` | 0 | include `<cblas.h>` and hand large problems to CBLAS |
| `LYGEMM_BLAS_THRESHOLD` | 32768 | `m*n*k` above which the runtime entry point uses BLAS |
| `LYGEMM_USE_JIT` | 1 on x86-64 Linux with AVX2 + FMA | generate a kernel per shape at run time |
| `LYGEMM_JIT_MAX_CODE` | 65536 | ceiling on the generated code per kernel, in bytes |

## Migration from the previous interface

`ly_dgemm` and the built library are gone; the project is header only.

```cpp
// before:  void ly_dgemm(int, int, int, const double*, const double*, double*)
//          built with -llygemm, and correct only for 4x4 and 2x4
// after:
#include <lygemm.hpp>
lygemm::gemm<4, 4, 4>(A, B, C);      // or lygemm::gemm(m, n, k, A, B, C)
```

The old kernels assumed `k % 4 == 0` for `m == 4` and hard coded `k == 4` for
`m == 2`, so they read past the end of the inputs for other shapes. The new
entry points handle any shape.

## License

MIT, see `LICENSE`.
