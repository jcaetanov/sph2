import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# ============================================================
# PARÂMETROS
# ============================================================
np.random.seed(42)

n_particles = 5000
frames_max = 100
dt = 0.03

mass = 0.003          # massa por partícula
h = 0.10              # smoothing length
G = 0.9
c = 50.0              # "velocidade da luz" reescalada
k = 0.08              # constante da EoS
gamma_eos = 5/3
eps = 0.03            # softening gravitacional

bh_mass = 0.0
bh_formed = False
alive = None

hawking_particles = []   # lista de [pos(3,), vel(3,), vida]

# câmera / projeção
cam_ax = np.deg2rad(25)   # rotação em x
cam_ay = np.deg2rad(35)   # rotação em y
perspective = False
cam_dist = 4.0

# ============================================================
# INICIALIZAÇÃO 3D
# ============================================================
def sample_sphere(n, radius=1.0):
    """
    Amostra aproximadamente uniforme dentro de uma esfera 3D.
    """
    u = np.random.rand(n)
    costheta = 2*np.random.rand(n) - 1
    phi = 2*np.pi*np.random.rand(n)

    r = radius * u**(1/3)
    sintheta = np.sqrt(1 - costheta**2)

    x = r * sintheta * np.cos(phi)
    y = r * sintheta * np.sin(phi)
    z = r * costheta
    return np.column_stack((x, y, z))

pos = sample_sphere(n_particles, radius=0.9)

# rotação inicial para formar disco/espiral no colapso
vel = np.zeros_like(pos)
omega = 0.45
vel[:, 0] = -omega * pos[:, 1]
vel[:, 1] =  omega * pos[:, 0]
vel += 0.03 * np.random.randn(*vel.shape)

alive = np.ones(n_particles, dtype=bool)

# ============================================================
# KERNEL SPH 3D
# ============================================================
def W_vec_3d(r, h):
    q = r / h
    sigma = 1.0 / (np.pi * h**3)
    W = np.zeros_like(r)

    mask1 = (q >= 0) & (q < 1)
    mask2 = (q >= 1) & (q < 2)

    W[mask1] = sigma * (1 - 1.5*q[mask1]**2 + 0.75*q[mask1]**3)
    W[mask2] = sigma * 0.25 * (2 - q[mask2])**3
    return W

def grad_W_vec_3d(r_vec, h):
    r_norm = np.linalg.norm(r_vec, axis=2)
    q = r_norm / h
    sigma = 1.0 / (np.pi * h**3)

    grad_scalar = np.zeros_like(r_norm)

    mask1 = (q >= 0) & (q < 1)
    mask2 = (q >= 1) & (q < 2)

    grad_scalar[mask1] = sigma * (-3*q[mask1] + 2.25*q[mask1]**2) / h
    grad_scalar[mask2] = -0.75 * sigma * (2 - q[mask2])**2 / h

    r_safe = r_norm.copy()
    r_safe[r_safe == 0] = 1e-12

    grad_vec = grad_scalar[:, :, None] * r_vec / r_safe[:, :, None]
    return grad_vec

# ============================================================
# FÍSICA 3D
# ============================================================
def schwarzschild_radius(M):
    return 2 * G * M / (c**2 + 1e-12)

def compute_density_3d(pos):
    r_vec = pos[:, None, :] - pos[None, :, :]
    r_norm = np.linalg.norm(r_vec, axis=2)
    rho = mass * W_vec_3d(r_norm, h).sum(axis=1)
    return rho

def compute_pressure_acc_3d(pos, rho):
    P = k * rho**gamma_eos

    r_vec = pos[:, None, :] - pos[None, :, :]
    gradW = grad_W_vec_3d(r_vec, h)

    rho_safe = rho + 1e-12
    P_term = (P[:, None] / rho_safe[:, None]**2 +
              P[None, :] / rho_safe[None, :]**2)

    acc_pressure = -mass * (P_term[:, :, None] * gradW).sum(axis=1)
    return acc_pressure

def compute_self_gravity_3d(pos):
    r_vec = pos[:, None, :] - pos[None, :, :]
    r2 = np.sum(r_vec**2, axis=2) + eps**2
    np.fill_diagonal(r2, np.inf)

    acc = -G * mass * (r_vec / (r2[:, :, None]**1.5)).sum(axis=1)
    return acc

def compute_bh_acc_3d(pos, M_bh):
    if M_bh <= 0:
        return np.zeros_like(pos)

    r = np.linalg.norm(pos, axis=1)
    rs = schwarzschild_radius(M_bh)

    r_safe = np.maximum(r, rs + 1e-3)
    a_mag = -G * M_bh / ((r_safe - rs)**2 + eps**2)

    r_hat = pos / r_safe[:, None]
    return a_mag[:, None] * r_hat

def effective_radius(pos_alive):
    if len(pos_alive) == 0:
        return 0.0
    r = np.linalg.norm(pos_alive, axis=1)
    return np.percentile(r, 85)

def compute_total_acc(pos, alive_mask, bh_mass):
    acc = np.zeros_like(pos)

    if np.sum(alive_mask) == 0:
        return acc, np.zeros(pos.shape[0])

    p = pos[alive_mask]

    rho = compute_density_3d(p)
    acc_pressure = compute_pressure_acc_3d(p, rho)
    acc_gravity = compute_self_gravity_3d(p)
    acc_bh = compute_bh_acc_3d(p, bh_mass)

    acc_alive = acc_pressure + acc_gravity + acc_bh
    acc[alive_mask] = acc_alive

    rho_all = np.zeros(pos.shape[0])
    rho_all[alive_mask] = rho
    return acc, rho_all

# ============================================================
# INTEGRADOR
# ============================================================
def rk4_step_3d(pos, vel, alive_mask, bh_mass, dt):
    a1, rho1 = compute_total_acc(pos, alive_mask, bh_mass)

    pos2 = pos.copy()
    vel2 = vel.copy()
    pos2[alive_mask] += 0.5 * dt * vel[alive_mask]
    vel2[alive_mask] += 0.5 * dt * a1[alive_mask]

    a2, _ = compute_total_acc(pos2, alive_mask, bh_mass)

    pos3 = pos.copy()
    vel3 = vel.copy()
    pos3[alive_mask] += 0.5 * dt * vel2[alive_mask]
    vel3[alive_mask] += 0.5 * dt * a2[alive_mask]

    a3, _ = compute_total_acc(pos3, alive_mask, bh_mass)

    pos4 = pos.copy()
    vel4 = vel.copy()
    pos4[alive_mask] += dt * vel3[alive_mask]
    vel4[alive_mask] += dt * a3[alive_mask]

    a4, _ = compute_total_acc(pos4, alive_mask, bh_mass)

    pos_next = pos.copy()
    vel_next = vel.copy()

    pos_next[alive_mask] += (dt/6.0) * (
        vel[alive_mask] + 2*vel2[alive_mask] + 2*vel3[alive_mask] + vel4[alive_mask]
    )
    vel_next[alive_mask] += (dt/6.0) * (
        a1[alive_mask] + 2*a2[alive_mask] + 2*a3[alive_mask] + a4[alive_mask]
    )

    return pos_next, vel_next, rho1

# ============================================================
# HAWKING 3D
# ============================================================
def random_unit_vectors(n):
    v = np.random.randn(n, 3)
    norm = np.linalg.norm(v, axis=1, keepdims=True) + 1e-12
    return v / norm

def spawn_hawking_particles(M_bh, n_new=2):
    global hawking_particles

    if M_bh <= 0:
        return

    rs = schwarzschild_radius(M_bh)
    T_visual = 1.0 / (M_bh + 1e-6)

    dirs = random_unit_vectors(n_new)
    for d in dirs:
        p = d * (rs + 0.01)
        speed = 0.05 + 0.25 * T_visual   # amplificado visualmente
        v = d * speed + 0.01*np.random.randn(3)
        life = np.random.randint(40, 90)
        hawking_particles.append([p, v, life])

def update_hawking(dt):
    global hawking_particles
    new_hp = []
    for p, v, life in hawking_particles:
        p = p + dt * v
        life -= 1
        if life > 0:
            new_hp.append([p, v, life])
    hawking_particles = new_hp

# ============================================================
# PROJEÇÃO 3D -> 2D
# ============================================================
def rotate_points(points, ax_angle, ay_angle):
    cx, sx = np.cos(ax_angle), np.sin(ax_angle)
    cy, sy = np.cos(ay_angle), np.sin(ay_angle)

    Rx = np.array([
        [1, 0, 0],
        [0, cx, -sx],
        [0, sx,  cx]
    ])

    Ry = np.array([
        [ cy, 0, sy],
        [  0, 1,  0],
        [-sy, 0, cy]
    ])

    R = Ry @ Rx
    return points @ R.T

def project_points(points):
    pr = rotate_points(points, cam_ax, cam_ay)
    x, y, z = pr[:, 0], pr[:, 1], pr[:, 2]

    if perspective:
        scale = cam_dist / (cam_dist - z + 1e-6)
        xp = x * scale
        yp = y * scale
    else:
        xp, yp = x, y

    return np.column_stack((xp, yp)), z

# ============================================================
# FIGURA
# ============================================================
fig, ax = plt.subplots(figsize=(7, 7), facecolor='black')
ax.set_facecolor('black')
ax.set_xlim(-2.0, 2.0)
ax.set_ylim(-2.0, 2.0)
ax.set_aspect('equal')
ax.axis('off')

matter_scatter = ax.scatter([], [], s=6, c=[], cmap='inferno', vmin=0, vmax=1, alpha=0.9)
hawk_scatter = ax.scatter([], [], s=8, c='cyan', alpha=0.75)
bh_circle = plt.Circle((0, 0), 0.0, color='black', zorder=20)
photon_ring = plt.Circle((0, 0), 0.0, color='orange', fill=False, lw=1.2, alpha=0.8, zorder=19)

ax.add_patch(photon_ring)
ax.add_patch(bh_circle)

step_text = ax.text(
    0.02, 0.97, '',
    transform=ax.transAxes,
    color='white',
    fontsize=10,
    ha='left',
    va='top',
    bbox=dict(facecolor='black', alpha=0.35, edgecolor='none')
)

# ============================================================
# UPDATE
# ============================================================
def update(frame):
    global pos, vel, alive, bh_mass, bh_formed

    if np.sum(alive) > 0:
        pos, vel, rho = rk4_step_3d(pos, vel, alive, bh_mass, dt)
    else:
        rho = np.zeros(pos.shape[0])

    # critério de formação do BH
    visible_mass = np.sum(alive) * mass
    total_mass = bh_mass + visible_mass

    if np.sum(alive) > 10:
        Reff = effective_radius(pos[alive])
    else:
        Reff = 0.0

    rs_total = schwarzschild_radius(total_mass)

    if (not bh_formed) and (visible_mass > 0) and (Reff < 1.2 * rs_total):
        bh_formed = True
        bh_mass = max(0.15 * total_mass, bh_mass)

    # absorção pelo horizonte
    if bh_formed and np.sum(alive) > 0:
        rs = schwarzschild_radius(max(bh_mass, 1e-9))
        r_alive = np.linalg.norm(pos, axis=1)
        absorbed = alive & (r_alive < rs)

        if np.any(absorbed):
            bh_mass += np.sum(absorbed) * mass
            alive[absorbed] = False

        spawn_hawking_particles(bh_mass, n_new=3)
        update_hawking(dt)

        bh_circle.set_radius(rs)
        photon_ring.set_radius(2.6 * rs)
    else:
        rs = 0.0
        bh_circle.set_radius(0.0)
        photon_ring.set_radius(0.0)

    # matéria viva
    if np.sum(alive) > 0:
        proj_alive, z_alive = project_points(pos[alive])

        rho_alive = rho[alive]
        if np.max(rho_alive) > np.min(rho_alive):
            rho_norm = (rho_alive - np.min(rho_alive)) / (np.max(rho_alive) - np.min(rho_alive))
        else:
            rho_norm = np.zeros_like(rho_alive)

        matter_scatter.set_offsets(proj_alive)
        matter_scatter.set_array(rho_norm)
    else:
        matter_scatter.set_offsets(np.empty((0, 2)))
        matter_scatter.set_array(np.array([]))

    # hawking
    if len(hawking_particles) > 0:
        hp_pos = np.array([item[0] for item in hawking_particles])
        hp_proj, _ = project_points(hp_pos)
        hawk_scatter.set_offsets(hp_proj)
    else:
        hawk_scatter.set_offsets(np.empty((0, 2)))


    step_text.set_text(
        f"Step: {frame}/{frames_max}\n"
        f"BH: {'sim' if bh_formed else 'não'}\n"
        f"M_BH: {bh_mass:.3f}\n"
        f"r_s: {rs:.4f}"

    )

    if frame % 10 == 0:
        print(
            f"Frame {frame}/{frames_max} | "
            f"BH={'sim' if bh_formed else 'não'} | "
            f"M_BH={bh_mass:.4f} | rs={rs:.4f}"
        )

    return matter_scatter, hawk_scatter, bh_circle, photon_ring, step_text

ani = FuncAnimation(fig, update, frames=frames_max, interval=40, blit=True)

writer = PillowWriter(fps=24)
ani.save("black_hole_birth.gif", writer=writer)

print("GIF gerado: black_hole_birth.gif")
