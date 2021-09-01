# vk_launcher

- First we need to set up Taichi for development by [Developer installation](https://taichi.readthedocs.io/en/stable/dev_install.html).

- Under the taichi repo, execute `TAICHI_CMAKE_ARGS="-DCMAKE_CXX_COMPILER=clang++ -DTI_WITH_VULKAN=ON" python setup.py develop` to build the library with vulkan.

- Check out `examples` for runnable examples. 
```
python examples/simulation/mpm88.py

```

- Let's use a simple example to test our launcher. 

```
import taichi as ti
 
ti.init(arch=ti.vulkan)
 
x = ti.field(ti.i32, shape=8)
 
@ti.kernel
def run():
  for i in x:
    x[i] = i + 1

run()
mod = ti.aot.Module(ti.vulkan)
mod.add_kernel(run)
mod.save("./dumped_data/", "test")

```
- You will get `test_metadata.tcb` `test_metadata.txt` and `test_run_c4_0_k0008_vk_t00.spv` under `/TAICHI_REPO_DIR/dumped_data/`. There are the spirv binary files and some metadata files we need for the Android launcher.

- Open the Android Launcher using Android Studio
- We need to sync the dumped files manually. First save the dumped data to `./assets/dumped_data/`
```
adb push --sync assets/dumped_data /storage/emulated/0/Android/data/com.example.{PROJECT_NAME}/files/

```
- Finally, we can run the program on the phone and get the correct data from the root buffer `[1 2 3 4 5 6 7 8]`.

- We can also use `mpm88.py` as an example.
```

# MPM-MLS in 88 lines of Taichi code, originally created by @yuanming-hu
import taichi as ti

ti.init(arch=ti.vulkan)

n_particles = 5
n_grid = 10
dx = 1 / n_grid
dt = 2e-4

p_rho = 1
p_vol = (dx * 0.5)**2
p_mass = p_vol * p_rho
# gravity = 9.8
bound = 3
E = 400

x = ti.Vector.field(2, float, n_particles)
v = ti.Vector.field(2, float, n_particles)
C = ti.Matrix.field(2, 2, float, n_particles)
J = ti.field(float, n_particles)

grid_v = ti.Vector.field(2, float, (n_grid, n_grid))
grid_m = ti.field(float, (n_grid, n_grid))
gravity = 9.8

@ti.kernel
def substep():
    for i, j in grid_m:
        grid_v[i, j] = [0, 0]
        grid_m[i, j] = 0
    for p in x:
        Xp = x[p] / dx
        base = int(Xp - 0.5)
        fx = Xp - base
        w = [0.5 * (1.5 - fx)**2, 0.75 - (fx - 1)**2, 0.5 * (fx - 0.5)**2]
        stress = -dt * 4 * E * p_vol * (J[p] - 1) / dx**2
        affine = ti.Matrix([[stress, 0], [0, stress]]) + p_mass * C[p]
        for i, j in ti.static(ti.ndrange(3, 3)):
            offset = ti.Vector([i, j])
            dpos = (offset - fx) * dx
            weight = w[i].x * w[j].y
            grid_v[base + offset] += weight * (p_mass * v[p] + affine @ dpos)
            grid_m[base + offset] += weight * p_mass
    for i, j in grid_m:
        if grid_m[i, j] > 0:
            grid_v[i, j] /= grid_m[i, j]
        grid_v[i, j].y -= dt * gravity
        if i < bound and grid_v[i, j].x < 0:
            grid_v[i, j].x = 0
        if i > n_grid - bound and grid_v[i, j].x > 0:
            grid_v[i, j].x = 0
        if j < bound and grid_v[i, j].y < 0:
            grid_v[i, j].y = 0
        if j > n_grid - bound and grid_v[i, j].y > 0:
            grid_v[i, j].y = 0
    for p in x:
        Xp = x[p] / dx
        base = int(Xp - 0.5)
        fx = Xp - base
        w = [0.5 * (1.5 - fx)**2, 0.75 - (fx - 1)**2, 0.5 * (fx - 0.5)**2]
        new_v = ti.Vector.zero(float, 2)
        new_C = ti.Matrix.zero(float, 2, 2)
        for i, j in ti.static(ti.ndrange(3, 3)):
            offset = ti.Vector([i, j])
            dpos = (offset - fx) * dx
            weight = w[i].x * w[j].y
            g_v = grid_v[base + offset]
            new_v += weight * g_v
            new_C += 4 * weight * g_v.outer_product(dpos) / dx**2
        v[p] = new_v
        x[p] += dt * v[p]
        J[p] *= 1 + dt * new_C.trace()
        C[p] = new_C

@ti.kernel
def init():
    x[0] = [0.29083378, 0.3435237]
    x[1] = [0.4307647, 0.57781583]
    x[2] = [0.417721, 0.28510436]
    x[3] = [0.5669118, 0.27478838]
    x[4] = [0.24649793, 0.24536438]
    for i in range(n_particles):
        v[i] = [0, -1]
        J[i] = 1

init()
substep()

m = ti.aot.Module(ti.vulkan)
m.add_kernel(init)
m.add_kernel(substep)
m.save('./dumped_data/', 'mpm88')
```
- Again, we will get `mpm88_metadata.tcb` `mpm88_metadata.txt` and `mpm88_substep_c4_0_k0016_vk_t00.spv` `mpm88_substep_c4_0_k0016_vk_t01.spv` `mpm88_substep_c4_0_k0016_vk_t02.spv` `mpm88_substep_c4_0_k0016_vk_t03.spv` `mpm88_init_c6_0_k0009_vk_t00.spv` `mpm88_init_c6_0_k0009_vk_t01.spv` under `/TAICHI_REPO_DIR/dumped_data/`. There are the spirv binary files and some metadata files we need for the Android launcher.
- Upload these files to our android phone
- Finally, we can run the program on the phone and get the correct data from the root buffer ` 
[[0.2909014  0.3436858 ]
 [0.4307647  0.5778154 ]
 [0.41768864 0.2851951 ]
 [0.56704193 0.27491662]
 [0.24659473 0.24539314]]`.
