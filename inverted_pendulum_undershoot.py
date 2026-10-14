"""
台車倒立擺：單一實數右半平面零點造成的步階下衝 (undershoot)

開迴路（力 → 台車位置）不穩定，但有一個實數 RHP 零點 s = +sqrt(g/L)。
狀態回授把四個極點都拉到左半平面，零點位置不變。
因此「位置命令 r → 台車位置 x」這個穩定、嚴格真分式的閉迴路
恰有 1 個（奇數個）實數 RHP 零點。依作業第 46 題，單位步階必有下衝。

本程式用非線性運動方程式積分，控制器是對線性化模型做極點配置。
角度約 15 度以內，線性設計與非線性軌跡幾乎重合。

參數
    M = 1.0 kg, m = 0.25 kg, L = 1.5 m, g = 9.81 m/s^2
    RHP 零點  z = +sqrt(g/L) = +2.557
    閉迴路極點  -3.5, -4.0, -4.8, -5.6
    閉迴路轉移函數（r → x，直流增益 = 1）
        T(s) = (-57.541 s^2 + 376.32) / (s^4 + 17.9 s^3 + 118.88 s^2 + 347.2 s + 376.32)
             = -57.541 (s - 2.557)(s + 2.557) / (穩定多項式)
"""

from vpython import *

# ---------- 物理參數 ----------
M = 1.0
m_bob = 0.25
rod_len = 1.5
grav = 9.81
r_cmd = 1.0  # 正向步階：目標台車位置 (m)

# u = -kx (x - r) - kv v - kth theta - kw omega
# 極點配置於 -3.5, -4, -4.8, -5.6。kx < 0：位置迴路的符號與一般負回授相反。
kx = -57.5412844036
kv = -53.0886850152
kth = -276.8944266053
kw = -106.4830275228

z_rhp = sqrt(grav / rod_len)

# ---------- 非線性運動方程式 ----------
def accelerations(theta, omega, force):
    """由目前角度、角速度與外力，解出台車加速度與角加速度。"""
    c = cos(theta)
    s = sin(theta)
    # (M+m) xdd + m L thdd cos = F + m L w^2 sin
    # xdd cos + L thdd         = g sin
    a11 = M + m_bob
    a12 = m_bob * rod_len * c
    a21 = c
    a22 = rod_len
    b1 = force + m_bob * rod_len * omega * omega * s
    b2 = grav * s
    det = a11 * a22 - a12 * a21
    xdd = (b1 * a22 - a12 * b2) / det
    thdd = (a11 * b2 - b1 * a21) / det
    return xdd, thdd


def derivatives(state, force):
    x, v, theta, omega = state
    xdd, thdd = accelerations(theta, omega, force)
    return [v, xdd, omega, thdd]


def control_force(state, reference):
    x, v, theta, omega = state
    return -kx * (x - reference) - kv * v - kth * theta - kw * omega


def rk4_step(state, reference, dt):
    def f(s):
        return derivatives(s, control_force(s, reference))

    k1 = f(state)
    s2 = [state[i] + 0.5 * dt * k1[i] for i in range(4)]
    k2 = f(s2)
    s3 = [state[i] + 0.5 * dt * k2[i] for i in range(4)]
    k3 = f(s3)
    s4 = [state[i] + dt * k3[i] for i in range(4)]
    k4 = f(s4)
    return [state[i] + (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(4)]


# ---------- 畫面 ----------
scene = canvas(
    title="台車倒立擺：右半平面零點與步階下衝",
    width=980,
    height=560,
    center=vec(0.30, 0.85, 0),
    background=vec(0.10, 0.11, 0.14),
)
scene.range = 2.15
scene.forward = vec(-0.15, -0.25, -1)

for i in range(-8, 12):
    shade = 0.22 if i % 2 == 0 else 0.38
    box(pos=vec(i * 0.25 + 0.125, -0.06, 0), size=vec(0.25, 0.03, 0.70), color=vec(shade, shade, shade))

rail = box(pos=vec(0.4, 0.02, 0), size=vec(5.2, 0.025, 0.08), color=vec(0.75, 0.62, 0.25))

for mark, name, col in ((0.0, "0", color.white), (r_cmd, "+1 m 目標", color.green)):
    box(pos=vec(mark, 0.18, 0), size=vec(0.015, 0.36, 0.015), color=col)
    label(pos=vec(mark, -0.22, 0), text=name, height=14, box=False, color=col, opacity=0)

cart = box(pos=vec(0, 0.16, 0), size=vec(0.36, 0.16, 0.28), color=vec(0.20, 0.55, 0.95))
wheel_y = 0.07
wheel_z = 0.16
wheel_r = 0.07
wheels = []
for sx in (-0.12, 0.12):
    for sz in (-wheel_z, wheel_z):
        wheels.append(cylinder(pos=vec(sx, wheel_y, sz), axis=vec(0, 0, 0.04 if sz > 0 else -0.04),
                               radius=wheel_r, color=vec(0.15, 0.15, 0.15)))

pivot_h = 0.26
rod = cylinder(pos=vec(0, pivot_h, 0), axis=vec(0, rod_len, 0), radius=0.018, color=vec(0.95, 0.85, 0.45))
bob = sphere(pos=vec(0, pivot_h + rod_len, 0), radius=0.075, color=vec(0.90, 0.25, 0.22))
cart_trail = curve(color=vec(0.35, 0.75, 1.0), radius=0.012)
bob_trail = curve(color=vec(1.0, 0.72, 0.25), radius=0.008)
force_arrow = arrow(pos=vec(0, 0.16, 0.22), axis=vec(-0.4, 0, 0), shaftwidth=0.03,
                    color=color.red, visible=False)
min_marker = box(pos=vec(0, 0.05, 0), size=vec(0.02, 0.10, 0.02), color=color.red, visible=False)

hud = label(pos=vec(-1.15, 1.85, 0), text="", height=16, box=False, align="left", color=color.white, opacity=0)

pos_graph = graph(
    title="台車位置 x(t)：正向步階命令下的下衝",
    xtitle="t (s)，步階發生在 t = 0",
    ytitle="x (m)",
    width=980,
    height=300,
    xmin=-0.4,
    xmax=4.5,
    ymin=-0.55,
    ymax=1.25,
    fast=False,
)
cmd_curve = gcurve(color=color.green, label="目標 r(t)", width=2)
zero_curve = gcurve(color=vec(0.65, 0.65, 0.65), label="x = 0", width=1)
x_curve = gcurve(color=vec(0.25, 0.70, 1.0), label="台車 x(t)", width=3)
min_dot = gdots(color=color.red, radius=6, label="下衝最低點")

ang_graph = graph(
    title="擺角 θ(t)（由鉛直向上量起，正值表示擺錘偏向 +x）",
    xtitle="t (s)",
    ytitle="θ (deg)",
    width=980,
    height=240,
    xmin=-0.4,
    xmax=4.5,
    ymin=-20,
    ymax=20,
    fast=False,
)
ang_curve = gcurve(color=vec(0.95, 0.55, 0.20), label="θ(t)", width=2)
ang_zero = gcurve(color=vec(0.65, 0.65, 0.65), width=1)

# ---------- 模擬狀態 ----------
PHYS_DT = 0.0005
SUBSTEPS = 8          # 每幀 0.004 s；rate(40) 約為 6 倍慢動作
SLOW_RATE = 40
state = [0.0, 0.0, 0.0, 0.0]
t_sim = -0.25         # 先停在原點，t = 0 才下步階
hold_frames = 0
saw_undershoot = False
finished = False
x_min = 0.0
t_at_min = 0.0


def caption_for(t, x, v, theta, force):
    if t < 0:
        return "初始靜止。即將在 t = 0 下達正向步階 r = +1 m。"
    if not saw_undershoot and force < 0 and x <= 0:
        return "步階剛下達：控制器先把台車推向 −x（錯誤方向）。擺錘因慣性留在右側，θ 變成正值。"
    if saw_undershoot and x < 0:
        return "已到達下衝谷底。擺子倒向目標後，台車改往 +x 去追它。"
    if x < r_cmd - 0.05:
        return "台車通過原點，朝 +1 m 前進，同時把擺拉回直立。"
    return "閉迴路穩定在目標。直流增益為 1，但途中必須先下衝。"


def apply_visual(state, force):
    x, v, theta, omega = state
    cart.pos.x = x
    for w, sx0 in zip(wheels, (-0.12, -0.12, 0.12, 0.12)):
        w.pos.x = x + sx0
    pivot = vec(x, pivot_h, 0)
    rod.pos = pivot
    rod.axis = vec(rod_len * sin(theta), rod_len * cos(theta), 0)
    bob.pos = pivot + rod.axis
    cart_trail.append(vec(x, 0.05, 0))
    bob_trail.append(vec(bob.pos.x, 0.05, 0))
    # 力箭頭：只顯示方向與相對大小，避免初始力把箭頭拉出畫面
    force_arrow.visible = abs(force) > 1.0
    force_arrow.pos = vec(x, 0.16, 0.22)
    shown = max(min(force, 40.0), -40.0) / 40.0 * 0.55
    force_arrow.axis = vec(shown, 0, 0)
    force_arrow.color = color.red if force < 0 else color.orange


def reset_sim():
    global state, t_sim, hold_frames, saw_undershoot, finished, x_min, t_at_min
    state = [0.0, 0.0, 0.0, 0.0]
    t_sim = -0.25
    hold_frames = 0
    saw_undershoot = False
    finished = False
    x_min = 0.0
    t_at_min = 0.0
    cart_trail.clear()
    bob_trail.clear()
    x_curve.data = []
    cmd_curve.data = []
    zero_curve.data = []
    ang_curve.data = []
    ang_zero.data = []
    min_dot.data = []
    min_marker.visible = False
    apply_visual(state, 0.0)


button(text="重播", bind=lambda _: reset_sim())
scene.append_to_caption(
    "\n青色軌跡是台車，橙色軌跡是擺錘在地面的投影。"
    "目標在右側 +1 m，但步階一下達台車會先往左。"
    "\n閉迴路零點：s = ±%.3f（只有正的那個在右半平面）。極點：−3.5, −4, −4.8, −5.6。\n"
    % z_rhp
)

reset_sim()
print("RHP zero z = +%.4f rad/s" % z_rhp)
print("Closed-loop step must undershoot: one real RHP zero, stable, strictly proper.")

while True:
    rate(SLOW_RATE)
    if hold_frames > 0:
        hold_frames -= 1
        continue
    if finished:
        continue

    reference = 0.0 if t_sim < 0.0 else r_cmd
    for _ in range(SUBSTEPS):
        prev_v = state[1]
        force = control_force(state, reference)
        state = rk4_step(state, reference, PHYS_DT)
        t_sim += PHYS_DT
        reference = 0.0 if t_sim < 0.0 else r_cmd
        if t_sim >= 0.0 and state[0] < x_min:
            x_min = state[0]
            t_at_min = t_sim
        if (not saw_undershoot) and t_sim > 0.05 and prev_v < 0.0 and state[1] >= 0.0 and state[0] < 0.0:
            saw_undershoot = True
            hold_frames = 50
            min_marker.visible = True
            min_marker.pos = vec(state[0], 0.08, 0)
            min_dot.plot(t_sim, state[0])

    x, v, theta, omega = state
    force = control_force(state, reference)
    apply_visual(state, force if t_sim >= 0.0 else 0.0)
    cmd_curve.plot(t_sim, reference)
    zero_curve.plot(t_sim, 0)
    x_curve.plot(t_sim, x)
    ang_zero.plot(t_sim, 0)
    ang_curve.plot(t_sim, theta * 180.0 / pi)

    hud.text = (
        "t = %5.2f s\n"
        "x = %+.3f m\n"
        "θ = %+.1f°\n"
        "F = %+.1f N\n"
        "x_min = %+.3f m\n\n%s"
        % (t_sim, x, theta * 180.0 / pi, force, x_min, caption_for(t_sim, x, v, theta, force))
    )

    if t_sim > 4.2 and abs(x - r_cmd) < 0.02 and abs(v) < 0.05 and abs(theta) < 0.05:
        finished = True
        hud.text = (
            "穩定於目標。\n"
            "x(∞) = +1 m，但 x 先降到 %+.3f m。\n"
            "這就是單一實數 RHP 零點的下衝。\n\n"
            "可按「重播」。" % x_min
        )
