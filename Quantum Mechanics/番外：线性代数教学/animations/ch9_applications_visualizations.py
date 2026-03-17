"""
第9章：应用专题 - Manim 可视化动画
===========================================

文件用途：为线性代数教学第9章制作4个可视化场景
输出格式：MP4 动画
运行命令：
  conda run -n teaching_env manim -ql --no_latex_cleanup ch9_applications_visualizations.py PortfolioEfficientFrontier
  conda run -n teaching_env manim -ql --no_latex_cleanup ch9_applications_visualizations.py TennisRacketEffect
  conda run -n teaching_env manim -ql --no_latex_cleanup ch9_applications_visualizations.py StressMohrCircle
  conda run -n teaching_env manim -ql --no_latex_cleanup ch9_applications_visualizations.py CovarianceEllipseAndPCA

作者：kyksj-1
"""

from manim import *
import numpy as np


# ============================================================
# 辅助函数
# ============================================================
def create_info_panel(title_text, title_color, items, position=UL):
    """
    创建信息面板（标题 + 若干行文字/公式）
    返回 VGroup
    """
    title = Text(title_text, font_size=26, color=title_color)
    group = VGroup(title)
    for item in items:
        group.add(item)
    group.arrange(DOWN, buff=0.2, aligned_edge=LEFT)
    box = SurroundingRectangle(group, color=title_color, buff=0.15,
                                fill_opacity=0.06, stroke_width=1.5)
    result = VGroup(group, box)
    result.to_corner(position).shift(DOWN * 0.3)
    return result


# ============================================================
# Scene 1: Markowitz 有效前沿 (Efficient Frontier)
# 在 sigma-mu 平面上展示：
#   - 多只资产的风险-收益散点
#   - 有效前沿（双曲线的上半部分）
#   - 全局最小方差组合 (GMV) 用星标标记
# ============================================================
class PortfolioEfficientFrontier(Scene):
    def construct(self):
        # ---------- 阶段1: 标题 ----------
        title = Text("Markowitz Efficient Frontier", font_size=32).to_edge(UP)
        subtitle = MathTex(
            r"\min_{\mathbf{w}} \; \mathbf{w}^T \Sigma \mathbf{w}"
            r"\quad \text{s.t.} \quad"
            r"\mathbf{w}^T \boldsymbol{\mu} = \mu_p,\;"
            r"\mathbf{w}^T \mathbf{1} = 1",
            font_size=26,
        ).next_to(title, DOWN, buff=0.3)

        self.play(Write(title), run_time=1)
        self.play(Write(subtitle), run_time=1.5)
        self.wait(1)
        self.play(
            FadeOut(subtitle),
            title.animate.scale(0.8).to_corner(UL),
            run_time=0.8,
        )

        # ---------- 阶段2: 坐标系 ----------
        axes = Axes(
            x_range=[0, 0.35, 0.05],
            y_range=[0, 0.20, 0.02],
            x_length=9,
            y_length=5.5,
            axis_config={
                "color": GREY_B,
                "stroke_width": 1.5,
                "include_numbers": False,
            },
            tips=True,
        ).shift(DOWN * 0.3 + RIGHT * 0.3)

        # 手动添加刻度标签（避免自动数字过于密集）
        x_ticks = VGroup()
        for val in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
            tick_label = Text(f"{val:.0%}", font_size=16, color=GREY_A)
            tick_label.next_to(axes.c2p(val, 0), DOWN, buff=0.15)
            x_ticks.add(tick_label)

        y_ticks = VGroup()
        for val in [0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18]:
            tick_label = Text(f"{val:.0%}", font_size=16, color=GREY_A)
            tick_label.next_to(axes.c2p(0, val), LEFT, buff=0.15)
            y_ticks.add(tick_label)

        x_axis_label = Text("Volatility (sigma)", font_size=22,
                            color=GREY_A).next_to(axes.x_axis, DOWN, buff=0.4)
        y_axis_label = Text("Expected Return (mu)", font_size=22,
                            color=GREY_A).rotate(PI / 2).next_to(
                                axes.y_axis, LEFT, buff=0.5)

        self.play(
            Create(axes), Write(x_axis_label), Write(y_axis_label),
            run_time=1.2,
        )
        self.play(
            LaggedStart(*[FadeIn(t) for t in x_ticks], lag_ratio=0.05),
            LaggedStart(*[FadeIn(t) for t in y_ticks], lag_ratio=0.05),
            run_time=0.8,
        )

        # ---------- 阶段3: 资产散点 ----------
        # 定义5只资产的 (sigma, mu)
        assets = [
            {"name": "Bonds", "sigma": 0.05, "mu": 0.03, "color": BLUE_C},
            {"name": "Large Cap", "sigma": 0.15, "mu": 0.09, "color": GREEN_C},
            {"name": "Mid Cap", "sigma": 0.20, "mu": 0.11, "color": TEAL_C},
            {"name": "Small Cap", "sigma": 0.25, "mu": 0.13, "color": ORANGE},
            {"name": "Emerging", "sigma": 0.30, "mu": 0.15, "color": RED_C},
        ]

        asset_dots = VGroup()
        asset_labels = VGroup()
        for a in assets:
            dot = Dot(axes.c2p(a["sigma"], a["mu"]),
                      color=a["color"], radius=0.08)
            label = Text(a["name"], font_size=16, color=a["color"]).next_to(
                dot, UR, buff=0.1)
            asset_dots.add(dot)
            asset_labels.add(label)

        self.play(
            LaggedStart(
                *[FadeIn(d, scale=1.5) for d in asset_dots],
                lag_ratio=0.2,
            ),
            run_time=1.5,
        )
        self.play(
            LaggedStart(*[FadeIn(l) for l in asset_labels], lag_ratio=0.15),
            run_time=1,
        )
        self.wait(1)

        # ---------- 阶段4: 有效前沿曲线 ----------
        # 用参数化的双曲线上半部分来近似有效前沿
        # 有效前沿: sigma(mu) 是 mu 的双曲函数
        # 参数: a, b, c 来自 sigma^2 = a*mu^2 - 2*b*mu + c
        # 选择合理参数使前沿通过资产区域附近
        a_coeff = 15.0
        b_coeff = 1.2
        c_coeff = 0.08

        # GMV 点: mu_gmv = b/a, sigma_gmv = sqrt(c - b^2/a)
        mu_gmv = b_coeff / a_coeff
        sigma_gmv_sq = c_coeff - b_coeff**2 / a_coeff
        sigma_gmv = np.sqrt(max(sigma_gmv_sq, 0.001))

        # 参数化有效前沿（上半部分: mu >= mu_gmv）
        frontier_pts = []
        for mu_val in np.linspace(mu_gmv, 0.19, 200):
            sigma_sq = a_coeff * mu_val**2 - 2 * b_coeff * mu_val + c_coeff
            if sigma_sq > 0:
                sigma_val = np.sqrt(sigma_sq)
                if 0 < sigma_val < 0.34:
                    frontier_pts.append(axes.c2p(sigma_val, mu_val))

        # 无效前沿（下半部分: mu < mu_gmv）
        inefficient_pts = []
        for mu_val in np.linspace(0.01, mu_gmv, 100):
            sigma_sq = a_coeff * mu_val**2 - 2 * b_coeff * mu_val + c_coeff
            if sigma_sq > 0:
                sigma_val = np.sqrt(sigma_sq)
                if 0 < sigma_val < 0.34:
                    inefficient_pts.append(axes.c2p(sigma_val, mu_val))

        # 绘制有效前沿
        if len(frontier_pts) > 2:
            frontier_curve = VMobject()
            frontier_curve.set_points_smoothly(frontier_pts)
            frontier_curve.set_stroke(YELLOW, width=3.5, opacity=1.0)

        # 绘制无效前沿（虚线，较浅）
        if len(inefficient_pts) > 2:
            inefficient_curve = DashedVMobject(
                VMobject().set_points_smoothly(inefficient_pts),
                num_dashes=30,
            )
            inefficient_curve.set_stroke(GREY, width=2, opacity=0.5)

        frontier_label = Text("Efficient Frontier", font_size=20,
                              color=YELLOW).next_to(
                                  axes.c2p(0.08, 0.17), RIGHT, buff=0.2)

        # 先画无效前沿（较低调），再画有效前沿
        self.play(Create(inefficient_curve), run_time=1.5)
        self.play(Create(frontier_curve), run_time=2)
        self.play(Write(frontier_label), run_time=0.8)
        self.wait(1)

        # ---------- 阶段5: GMV 点标记 ----------
        gmv_screen = axes.c2p(sigma_gmv, mu_gmv)
        # 星标用一个稍大的点 + 外圈表示
        gmv_star = Star(n=5, outer_radius=0.18, inner_radius=0.09,
                        color=RED, fill_opacity=1.0).move_to(gmv_screen)
        gmv_label = Text("GMV", font_size=20, color=RED).next_to(
            gmv_star, DL, buff=0.15)
        gmv_detail = MathTex(
            r"\sigma_{\min} = " + f"{sigma_gmv:.1%}",
            font_size=22, color=RED,
        ).next_to(gmv_label, DOWN, buff=0.1)

        # 从 GMV 画虚线到两轴
        gmv_h_line = DashedLine(
            axes.c2p(0, mu_gmv), gmv_screen,
            color=RED, stroke_width=1.5,
        )
        gmv_v_line = DashedLine(
            gmv_screen, axes.c2p(sigma_gmv, 0),
            color=RED, stroke_width=1.5,
        )

        self.play(
            FadeIn(gmv_star, scale=2),
            Create(gmv_h_line), Create(gmv_v_line),
            run_time=1.2,
        )
        self.play(Write(gmv_label), Write(gmv_detail), run_time=0.8)
        self.wait(1)

        # ---------- 阶段6: 可行域着色（可选：淡色填充） ----------
        # 用半透明填充前沿内侧区域
        # 生成包围区域的点
        fill_pts = list(frontier_pts) + list(reversed(inefficient_pts))
        if len(fill_pts) > 4:
            fill_region = VMobject()
            fill_region.set_points_smoothly(fill_pts)
            fill_region.set_fill(YELLOW, opacity=0.08)
            fill_region.set_stroke(width=0)
            self.play(FadeIn(fill_region), run_time=1)

        # ---------- 阶段7: 总结 ----------
        summary = Text(
            "Covariance matrix determines portfolio risk structure",
            font_size=22, color=GREEN_A,
        ).to_edge(DOWN)
        self.play(FadeIn(summary), run_time=1)
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 2: 网球拍效应 (Tennis Racket Effect / Intermediate Axis Theorem)
# 展示中间轴旋转的不稳定性：
#   - 三个主惯量轴用不同颜色
#   - 在三个投影平面上画相轨迹
#   - 轴1和轴3附近的轨迹闭合（稳定），轴2附近发散（不稳定）
# ============================================================
class TennisRacketEffect(Scene):
    def construct(self):
        # ---------- 阶段1: 标题和物理背景 ----------
        title = Text("Tennis Racket Effect (Intermediate Axis Theorem)",
                      font_size=28).to_edge(UP)
        # 欧拉方程
        euler_eq = MathTex(
            r"I_1 \dot{\omega}_1 = (I_2 - I_3)\omega_2 \omega_3",
            font_size=28,
        ).shift(UP * 2)
        euler_eq2 = MathTex(
            r"I_2 \dot{\omega}_2 = (I_3 - I_1)\omega_3 \omega_1",
            font_size=28,
        ).next_to(euler_eq, DOWN, buff=0.2)
        euler_eq3 = MathTex(
            r"I_3 \dot{\omega}_3 = (I_1 - I_2)\omega_1 \omega_2",
            font_size=28,
        ).next_to(euler_eq2, DOWN, buff=0.2)

        self.play(Write(title), run_time=1)
        self.play(
            Write(euler_eq), Write(euler_eq2), Write(euler_eq3),
            run_time=2,
        )
        self.wait(1.5)

        # 惯量标注
        inertia_label = MathTex(
            r"I_1 < I_2 < I_3",
            font_size=32, color=YELLOW,
        ).next_to(euler_eq3, DOWN, buff=0.4)
        inertia_box = SurroundingRectangle(inertia_label, color=YELLOW, buff=0.12)
        self.play(Write(inertia_label), Create(inertia_box), run_time=1)
        self.wait(1)

        # 清理
        self.play(
            FadeOut(euler_eq), FadeOut(euler_eq2), FadeOut(euler_eq3),
            FadeOut(inertia_label), FadeOut(inertia_box),
            title.animate.scale(0.75).to_corner(UL),
            run_time=0.8,
        )

        # ---------- 阶段2: 惯量值和物理设置 ----------
        I1, I2, I3 = 1.0, 2.0, 3.0  # I1 < I2 < I3

        # 定义欧拉方程的右端函数
        def euler_rhs(omega, I1, I2, I3):
            w1, w2, w3 = omega
            dw1 = (I2 - I3) * w2 * w3 / I1
            dw2 = (I3 - I1) * w3 * w1 / I2
            dw3 = (I1 - I2) * w1 * w2 / I3
            return np.array([dw1, dw2, dw3])

        # RK4 积分
        def integrate_euler(omega0, I1, I2, I3, dt=0.005, steps=4000):
            traj = [omega0.copy()]
            omega = omega0.copy()
            for _ in range(steps):
                k1 = euler_rhs(omega, I1, I2, I3)
                k2 = euler_rhs(omega + 0.5 * dt * k1, I1, I2, I3)
                k3 = euler_rhs(omega + 0.5 * dt * k2, I1, I2, I3)
                k4 = euler_rhs(omega + dt * k3, I1, I2, I3)
                omega = omega + dt / 6.0 * (k1 + 2*k2 + 2*k3 + k4)
                traj.append(omega.copy())
            return np.array(traj)

        # ---------- 阶段3: 三个子图并排 ----------
        # 在同一画面上展示三个相图
        # 左: omega2-omega3 平面 (绕轴1旋转)
        # 中: omega1-omega3 平面 (绕轴2旋转)
        # 右: omega1-omega2 平面 (绕轴3旋转)

        sub_width = 3.5
        sub_height = 3.5
        spacing = 0.6

        # 三个坐标系
        axes1 = Axes(
            x_range=[-1.5, 1.5, 0.5],
            y_range=[-1.5, 1.5, 0.5],
            x_length=sub_width,
            y_length=sub_height,
            axis_config={"color": GREY_B, "stroke_width": 1.2,
                         "include_numbers": False},
            tips=False,
        ).shift(LEFT * (sub_width + spacing) + DOWN * 0.3)

        axes2 = Axes(
            x_range=[-1.5, 1.5, 0.5],
            y_range=[-1.5, 1.5, 0.5],
            x_length=sub_width,
            y_length=sub_height,
            axis_config={"color": GREY_B, "stroke_width": 1.2,
                         "include_numbers": False},
            tips=False,
        ).shift(DOWN * 0.3)

        axes3 = Axes(
            x_range=[-1.5, 1.5, 0.5],
            y_range=[-1.5, 1.5, 0.5],
            x_length=sub_width,
            y_length=sub_height,
            axis_config={"color": GREY_B, "stroke_width": 1.2,
                         "include_numbers": False},
            tips=False,
        ).shift(RIGHT * (sub_width + spacing) + DOWN * 0.3)

        # 子图标题
        label1 = VGroup(
            Text("Axis 1 (stable)", font_size=18, color=BLUE_C),
            MathTex(r"I_1 = 1", font_size=20, color=BLUE_C),
        ).arrange(DOWN, buff=0.1).next_to(axes1, UP, buff=0.15)

        label2 = VGroup(
            Text("Axis 2 (UNSTABLE)", font_size=18, color=RED_C),
            MathTex(r"I_2 = 2", font_size=20, color=RED_C),
        ).arrange(DOWN, buff=0.1).next_to(axes2, UP, buff=0.15)

        label3 = VGroup(
            Text("Axis 3 (stable)", font_size=18, color=GREEN_C),
            MathTex(r"I_3 = 3", font_size=20, color=GREEN_C),
        ).arrange(DOWN, buff=0.1).next_to(axes3, UP, buff=0.15)

        # 坐标轴标签
        ax1_xlabel = MathTex(r"\omega_2", font_size=18).next_to(
            axes1.x_axis.get_end(), DR, buff=0.05)
        ax1_ylabel = MathTex(r"\omega_3", font_size=18).next_to(
            axes1.y_axis.get_end(), UL, buff=0.05)

        ax2_xlabel = MathTex(r"\omega_1", font_size=18).next_to(
            axes2.x_axis.get_end(), DR, buff=0.05)
        ax2_ylabel = MathTex(r"\omega_3", font_size=18).next_to(
            axes2.y_axis.get_end(), UL, buff=0.05)

        ax3_xlabel = MathTex(r"\omega_1", font_size=18).next_to(
            axes3.x_axis.get_end(), DR, buff=0.05)
        ax3_ylabel = MathTex(r"\omega_2", font_size=18).next_to(
            axes3.y_axis.get_end(), UL, buff=0.05)

        self.play(
            Create(axes1), Create(axes2), Create(axes3),
            Write(label1), Write(label2), Write(label3),
            Write(ax1_xlabel), Write(ax1_ylabel),
            Write(ax2_xlabel), Write(ax2_ylabel),
            Write(ax3_xlabel), Write(ax3_ylabel),
            run_time=1.5,
        )

        # ---------- 阶段4: 计算并绘制轨迹 ----------
        # 绕轴1旋转的扰动（稳定）
        axis1_trajs = []
        perturb_magnitudes = [0.1, 0.2, 0.35, 0.5]
        for eps in perturb_magnitudes:
            # 主分量在 omega1 方向, 小扰动在 omega2, omega3
            omega0 = np.array([1.0, eps, eps * 0.5])
            traj = integrate_euler(omega0, I1, I2, I3, dt=0.005, steps=3000)
            axis1_trajs.append(traj)

        # 绕轴2旋转的扰动（不稳定）
        axis2_trajs = []
        for eps in [0.05, 0.1, 0.15, 0.2]:
            omega0 = np.array([eps, 1.0, eps * 0.5])
            traj = integrate_euler(omega0, I1, I2, I3, dt=0.005, steps=3000)
            axis2_trajs.append(traj)

        # 绕轴3旋转的扰动（稳定）
        axis3_trajs = []
        for eps in perturb_magnitudes:
            omega0 = np.array([eps, eps * 0.5, 1.0])
            traj = integrate_euler(omega0, I1, I2, I3, dt=0.005, steps=3000)
            axis3_trajs.append(traj)

        # 绘制轨迹的辅助函数
        def draw_trajectory(ax, traj, ix, iy, color, opacity=0.8):
            """在给定的axes上，用traj的第ix和iy分量画轨迹"""
            pts = []
            for row in traj:
                x_val, y_val = row[ix], row[iy]
                if abs(x_val) < 1.45 and abs(y_val) < 1.45:
                    pts.append(ax.c2p(x_val, y_val))
            if len(pts) > 5:
                curve = VMobject()
                curve.set_points_smoothly(pts)
                curve.set_stroke(color, width=1.8, opacity=opacity)
                return curve
            return None

        # 轴1: 在 (omega2, omega3) 平面上的轨迹 -> 闭合小圈
        axis1_curves = VGroup()
        colors1 = color_gradient([BLUE_A, BLUE_D], len(axis1_trajs))
        for traj, col in zip(axis1_trajs, colors1):
            curve = draw_trajectory(axes1, traj, 1, 2, col)
            if curve is not None:
                axis1_curves.add(curve)

        # 轴2: 在 (omega1, omega3) 平面上的轨迹 -> 发散
        axis2_curves = VGroup()
        colors2 = color_gradient([RED_A, RED_D], len(axis2_trajs))
        for traj, col in zip(axis2_trajs, colors2):
            curve = draw_trajectory(axes2, traj, 0, 2, col)
            if curve is not None:
                axis2_curves.add(curve)

        # 轴3: 在 (omega1, omega2) 平面上的轨迹 -> 闭合小圈
        axis3_curves = VGroup()
        colors3 = color_gradient([GREEN_A, GREEN_D], len(axis3_trajs))
        for traj, col in zip(axis3_trajs, colors3):
            curve = draw_trajectory(axes3, traj, 0, 1, col)
            if curve is not None:
                axis3_curves.add(curve)

        # 标记原点 (平衡点)
        eq_dot1 = Dot(axes1.c2p(0, 0), color=BLUE_C, radius=0.05)
        eq_dot2 = Dot(axes2.c2p(0, 0), color=RED_C, radius=0.05)
        eq_dot3 = Dot(axes3.c2p(0, 0), color=GREEN_C, radius=0.05)

        # 动画: 依次显示三个子图的轨迹
        self.play(
            FadeIn(eq_dot1), FadeIn(eq_dot2), FadeIn(eq_dot3),
            run_time=0.5,
        )

        # 轴1 的轨迹 (稳定)
        self.play(
            LaggedStart(*[Create(c) for c in axis1_curves], lag_ratio=0.2),
            run_time=2,
        )
        stable1_text = Text("Closed orbits", font_size=16,
                            color=BLUE_C).next_to(axes1, DOWN, buff=0.15)
        self.play(FadeIn(stable1_text), run_time=0.5)

        # 轴2 的轨迹 (不稳定)
        self.play(
            LaggedStart(*[Create(c) for c in axis2_curves], lag_ratio=0.2),
            run_time=2,
        )
        unstable_text = Text("Trajectories diverge!", font_size=16,
                             color=RED_C).next_to(axes2, DOWN, buff=0.15)
        self.play(FadeIn(unstable_text), run_time=0.5)

        # 轴3 的轨迹 (稳定)
        self.play(
            LaggedStart(*[Create(c) for c in axis3_curves], lag_ratio=0.2),
            run_time=2,
        )
        stable3_text = Text("Closed orbits", font_size=16,
                            color=GREEN_C).next_to(axes3, DOWN, buff=0.15)
        self.play(FadeIn(stable3_text), run_time=0.5)
        self.wait(1)

        # ---------- 阶段5: 总结 ----------
        summary = Text(
            "Intermediate axis rotation is unstable (eigenvalue analysis of linearized Euler equations)",
            font_size=20, color=YELLOW,
        ).to_edge(DOWN)
        self.play(FadeIn(summary), run_time=1)
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 3: Mohr 圆 (Mohr's Circle)
# 给定三个主应力 sigma_1 > sigma_2 > sigma_3
# 画三个嵌套 Mohr 圆
# 标注主应力点和最大剪应力
# ============================================================
class StressMohrCircle(Scene):
    def construct(self):
        # ---------- 阶段1: 标题和背景知识 ----------
        title = Text("Mohr's Circle for 3D Stress State", font_size=30).to_edge(UP)

        # 应力张量的谱分解
        stress_eq = MathTex(
            r"\boldsymbol{\sigma} = Q \,"
            r"\text{diag}(\sigma_1, \sigma_2, \sigma_3)"
            r"\, Q^T",
            font_size=30,
        ).next_to(title, DOWN, buff=0.3)

        self.play(Write(title), run_time=1)
        self.play(Write(stress_eq), run_time=1.2)
        self.wait(1)
        self.play(
            FadeOut(stress_eq),
            title.animate.scale(0.8).to_corner(UL),
            run_time=0.8,
        )

        # ---------- 阶段2: 定义主应力和坐标系 ----------
        # 主应力值
        s1, s2, s3 = 120.0, 50.0, -30.0  # sigma_1 > sigma_2 > sigma_3 (MPa)

        # 三个 Mohr 圆的参数 (圆心, 半径)
        # 圆1: sigma_1 和 sigma_3 -> 最大圆
        c13 = (s1 + s3) / 2.0
        r13 = (s1 - s3) / 2.0
        # 圆2: sigma_1 和 sigma_2
        c12 = (s1 + s2) / 2.0
        r12 = (s1 - s2) / 2.0
        # 圆3: sigma_2 和 sigma_3
        c23 = (s2 + s3) / 2.0
        r23 = (s2 - s3) / 2.0

        # 坐标系范围
        x_min = s3 - 20
        x_max = s1 + 20
        y_max = r13 + 15

        axes = Axes(
            x_range=[x_min, x_max, 20],
            y_range=[-y_max, y_max, 20],
            x_length=10,
            y_length=6,
            axis_config={
                "color": GREY_B,
                "stroke_width": 1.5,
                "include_numbers": False,
            },
            tips=True,
        ).shift(DOWN * 0.2)

        # 坐标轴标签
        x_label = Text("Normal stress (sigma_n)", font_size=20,
                        color=GREY_A).next_to(axes.x_axis, DOWN, buff=0.3)
        y_label = Text("Shear stress (tau_n)", font_size=20,
                        color=GREY_A).rotate(PI / 2).next_to(
                            axes.y_axis, LEFT, buff=0.4)

        # 刻度标注
        tick_labels = VGroup()
        for val in [-20, 0, 20, 40, 60, 80, 100, 120]:
            if x_min <= val <= x_max:
                lbl = Text(str(int(val)), font_size=14, color=GREY_A)
                lbl.next_to(axes.c2p(val, 0), DOWN, buff=0.1)
                tick_labels.add(lbl)

        self.play(
            Create(axes), Write(x_label), Write(y_label),
            run_time=1.2,
        )
        self.play(
            LaggedStart(*[FadeIn(t) for t in tick_labels], lag_ratio=0.03),
            run_time=0.5,
        )

        # ---------- 阶段3: 主应力点 ----------
        # 主应力在实轴上（剪应力为0）
        dot_s1 = Dot(axes.c2p(s1, 0), color=RED, radius=0.08)
        dot_s2 = Dot(axes.c2p(s2, 0), color=GREEN, radius=0.08)
        dot_s3 = Dot(axes.c2p(s3, 0), color=BLUE, radius=0.08)

        lbl_s1 = MathTex(r"\sigma_1 = 120", font_size=20,
                          color=RED).next_to(dot_s1, UP, buff=0.15)
        lbl_s2 = MathTex(r"\sigma_2 = 50", font_size=20,
                          color=GREEN).next_to(dot_s2, DOWN, buff=0.15)
        lbl_s3 = MathTex(r"\sigma_3 = -30", font_size=20,
                          color=BLUE).next_to(dot_s3, UP, buff=0.15)

        self.play(
            FadeIn(dot_s1, scale=1.5), FadeIn(dot_s2, scale=1.5),
            FadeIn(dot_s3, scale=1.5),
            Write(lbl_s1), Write(lbl_s2), Write(lbl_s3),
            run_time=1.5,
        )
        self.wait(0.5)

        # ---------- 阶段4: 绘制三个 Mohr 圆 ----------
        # 辅助函数：在 axes 上创建圆
        def make_mohr_circle(center_x, radius, color, opacity=0.7):
            """在 sigma_n - tau_n 平面上画一个 Mohr 圆"""
            pts = []
            for t in np.linspace(0, 2 * np.pi, 200):
                sx = center_x + radius * np.cos(t)
                ty = radius * np.sin(t)
                pts.append(axes.c2p(sx, ty))
            pts.append(pts[0])  # 闭合
            circle = VMobject()
            circle.set_points_smoothly(pts)
            circle.set_stroke(color, width=2.5, opacity=opacity)
            return circle

        # 最大圆 (sigma_1 与 sigma_3)
        circle_13 = make_mohr_circle(c13, r13, YELLOW, opacity=0.9)
        circle_label_13 = MathTex(
            r"\sigma_1\text{-}\sigma_3", font_size=18, color=YELLOW,
        ).move_to(axes.c2p(c13, r13 + 8))

        # 圆 (sigma_1 与 sigma_2)
        circle_12 = make_mohr_circle(c12, r12, RED_C, opacity=0.8)
        circle_label_12 = MathTex(
            r"\sigma_1\text{-}\sigma_2", font_size=18, color=RED_C,
        ).move_to(axes.c2p(c12, r12 + 8))

        # 圆 (sigma_2 与 sigma_3)
        circle_23 = make_mohr_circle(c23, r23, BLUE_C, opacity=0.8)
        circle_label_23 = MathTex(
            r"\sigma_2\text{-}\sigma_3", font_size=18, color=BLUE_C,
        ).move_to(axes.c2p(c23, r23 + 8))

        # 动画：先画最大圆，再两个小圆
        self.play(Create(circle_13), Write(circle_label_13), run_time=2)
        self.wait(0.5)
        self.play(
            Create(circle_12), Write(circle_label_12),
            Create(circle_23), Write(circle_label_23),
            run_time=2,
        )
        self.wait(0.5)

        # ---------- 阶段5: 标注最大剪应力 ----------
        tau_max = r13  # 最大剪应力 = (sigma_1 - sigma_3) / 2
        tau_max_dot = Dot(axes.c2p(c13, tau_max), color=YELLOW, radius=0.08)
        tau_max_label = MathTex(
            r"\tau_{\max} = \frac{\sigma_1 - \sigma_3}{2} = " + f"{tau_max:.0f}",
            font_size=22, color=YELLOW,
        ).next_to(tau_max_dot, UR, buff=0.15)

        # 虚线从最大剪应力点到两轴
        dashed_h = DashedLine(
            axes.c2p(0, tau_max), axes.c2p(c13, tau_max),
            color=YELLOW, stroke_width=1.5,
        )
        dashed_v = DashedLine(
            axes.c2p(c13, 0), axes.c2p(c13, tau_max),
            color=YELLOW, stroke_width=1.5,
        )

        self.play(
            FadeIn(tau_max_dot, scale=1.5),
            Create(dashed_h), Create(dashed_v),
            run_time=1,
        )
        self.play(Write(tau_max_label), run_time=0.8)
        self.wait(1)

        # ---------- 阶段6: 可行应力状态区域着色 ----------
        # 可行区域 = 最大圆内 且 在两个小圆外 的区域（用文字说明即可）
        feasible_text = Text(
            "All stress states (sigma_n, tau_n) lie in the shaded region",
            font_size=18, color=GREY_A,
        ).to_edge(DOWN)
        self.play(FadeIn(feasible_text), run_time=0.8)
        self.wait(1)

        # ---------- 阶段7: 总结 ----------
        summary_box_content = VGroup(
            MathTex(r"\sigma_1 > \sigma_2 > \sigma_3", font_size=24),
            MathTex(
                r"\tau_{\max} = \frac{\sigma_1 - \sigma_3}{2}",
                font_size=24, color=YELLOW,
            ),
        ).arrange(DOWN, buff=0.2).to_corner(DR).shift(UP * 0.5 + LEFT * 0.3)
        summary_rect = SurroundingRectangle(summary_box_content, color=WHITE,
                                             buff=0.15, fill_opacity=0.05)

        self.play(
            FadeIn(summary_box_content), Create(summary_rect),
            run_time=1,
        )
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 4: 协方差椭圆与 PCA (Covariance Ellipse and PCA)
# 展示：
#   - 有相关性的 2D 散点数据
#   - 协方差椭圆
#   - 两个主成分方向（特征向量箭头）
#   - 长箭头 = PC1, 短箭头 = PC2
# ============================================================
class CovarianceEllipseAndPCA(Scene):
    def construct(self):
        # ---------- 阶段1: 标题 ----------
        title = Text("Covariance Ellipse and PCA", font_size=32).to_edge(UP)
        pca_eq = MathTex(
            r"\Sigma = Q \Lambda Q^T"
            r"\quad \Longrightarrow \quad"
            r"\text{PC directions} = \text{eigenvectors of } \Sigma",
            font_size=26,
        ).next_to(title, DOWN, buff=0.3)

        self.play(Write(title), run_time=1)
        self.play(Write(pca_eq), run_time=1.5)
        self.wait(1)
        self.play(
            FadeOut(pca_eq),
            title.animate.scale(0.8).to_corner(UL),
            run_time=0.8,
        )

        # ---------- 阶段2: 生成数据和坐标系 ----------
        # 协方差矩阵
        cov_matrix = np.array([[2.5, 1.8], [1.8, 1.8]])
        mean = np.array([0, 0])

        # 生成数据
        np.random.seed(42)
        data = np.random.multivariate_normal(mean, cov_matrix, 200)

        # 特征分解
        evals, evecs = np.linalg.eigh(cov_matrix)
        # evals 升序, 反转为降序 (PC1 对应最大特征值)
        idx = np.argsort(evals)[::-1]
        evals = evals[idx]
        evecs = evecs[:, idx]

        # 坐标系
        plane = NumberPlane(
            x_range=[-5, 5, 1],
            y_range=[-4, 4, 1],
            x_length=9,
            y_length=6,
            background_line_style={
                "stroke_color": BLUE_E,
                "stroke_width": 0.6,
                "stroke_opacity": 0.2,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1.5},
        ).shift(DOWN * 0.2)

        x_lbl = plane.get_x_axis_label(MathTex("x_1", font_size=24))
        y_lbl = plane.get_y_axis_label(MathTex("x_2", font_size=24))

        self.play(Create(plane), Write(x_lbl), Write(y_lbl), run_time=1)

        # ---------- 阶段3: 散点数据 ----------
        dots = VGroup()
        for pt in data:
            if abs(pt[0]) < 4.8 and abs(pt[1]) < 3.8:
                d = Dot(plane.c2p(pt[0], pt[1]),
                        color=BLUE_C, radius=0.035, fill_opacity=0.6)
                dots.add(d)

        self.play(
            LaggedStart(*[FadeIn(d, scale=0.5) for d in dots], lag_ratio=0.005),
            run_time=2,
        )
        self.wait(0.5)

        # 信息面板：协方差矩阵
        cov_label = MathTex(
            r"\Sigma = \begin{pmatrix} 2.5 & 1.8 \\ 1.8 & 1.8 \end{pmatrix}",
            font_size=26,
        ).to_corner(UR).shift(DOWN * 0.5 + LEFT * 0.2)
        cov_box = SurroundingRectangle(cov_label, color=TEAL, buff=0.12,
                                        fill_opacity=0.06)

        self.play(Write(cov_label), Create(cov_box), run_time=1)
        self.wait(0.5)

        # ---------- 阶段4: 协方差椭圆 ----------
        # 椭圆: x^T Sigma^{-1} x = c, 取 c 使得覆盖约 95% 数据 (c = 5.991 for chi2(2))
        # 参数化: x = sqrt(c) * Q * diag(sqrt(lam)) * [cos(t), sin(t)]
        chi2_val = 5.991  # 95% 置信度
        ellipse_pts = []
        for t in np.linspace(0, 2 * np.pi, 200):
            u = np.array([np.cos(t), np.sin(t)])
            # 椭圆上的点
            pt = np.sqrt(chi2_val) * evecs @ (np.sqrt(evals) * u)
            if abs(pt[0]) < 4.8 and abs(pt[1]) < 3.8:
                ellipse_pts.append(plane.c2p(pt[0], pt[1]))
        ellipse_pts.append(ellipse_pts[0])  # 闭合

        ellipse = VMobject()
        ellipse.set_points_smoothly(ellipse_pts)
        ellipse.set_stroke(YELLOW, width=3, opacity=0.9)

        ellipse_label = Text("95% Covariance Ellipse", font_size=18,
                             color=YELLOW).next_to(
                                 plane.c2p(evecs[0, 0] * np.sqrt(evals[0]) * 2,
                                           evecs[1, 0] * np.sqrt(evals[0]) * 2),
                                 UR, buff=0.2)

        self.play(Create(ellipse), run_time=2)
        self.play(Write(ellipse_label), run_time=0.8)
        self.wait(1)

        # ---------- 阶段5: 主成分方向 ----------
        # PC1: 最大特征值方向（长箭头）
        pc1_dir = evecs[:, 0]
        pc1_len = np.sqrt(evals[0]) * 2.0  # 长度按标准差缩放
        pc1_arrow = Arrow(
            plane.c2p(0, 0),
            plane.c2p(*(pc1_dir * pc1_len)),
            buff=0, stroke_width=5, color=RED,
            max_tip_length_to_length_ratio=0.1,
        )
        pc1_label = MathTex(
            r"\text{PC1}: \lambda_1 = " + f"{evals[0]:.2f}",
            font_size=22, color=RED,
        ).next_to(pc1_arrow.get_end(), UR, buff=0.15)

        # PC2: 最小特征值方向（短箭头）
        pc2_dir = evecs[:, 1]
        pc2_len = np.sqrt(evals[1]) * 2.0
        pc2_arrow = Arrow(
            plane.c2p(0, 0),
            plane.c2p(*(pc2_dir * pc2_len)),
            buff=0, stroke_width=5, color=GREEN,
            max_tip_length_to_length_ratio=0.15,
        )
        pc2_label = MathTex(
            r"\text{PC2}: \lambda_2 = " + f"{evals[1]:.2f}",
            font_size=22, color=GREEN,
        ).next_to(pc2_arrow.get_end(), DL, buff=0.15)

        # 反方向的虚线（显示完整的方向轴）
        pc1_dashed = DashedLine(
            plane.c2p(0, 0),
            plane.c2p(*(- pc1_dir * pc1_len)),
            color=RED, stroke_width=2, stroke_opacity=0.4,
        )
        pc2_dashed = DashedLine(
            plane.c2p(0, 0),
            plane.c2p(*(- pc2_dir * pc2_len)),
            color=GREEN, stroke_width=2, stroke_opacity=0.4,
        )

        self.play(
            GrowArrow(pc1_arrow), Create(pc1_dashed),
            run_time=1.5,
        )
        self.play(Write(pc1_label), run_time=0.8)

        self.play(
            GrowArrow(pc2_arrow), Create(pc2_dashed),
            run_time=1.5,
        )
        self.play(Write(pc2_label), run_time=0.8)
        self.wait(1)

        # ---------- 阶段6: 方差解释比 ----------
        total_var = evals.sum()
        var_ratio_1 = evals[0] / total_var * 100
        var_ratio_2 = evals[1] / total_var * 100

        var_info = VGroup(
            Text("Variance explained:", font_size=20, color=GREY_A),
            MathTex(
                r"\text{PC1}: " + f"{var_ratio_1:.1f}" + r"\%",
                font_size=22, color=RED,
            ),
            MathTex(
                r"\text{PC2}: " + f"{var_ratio_2:.1f}" + r"\%",
                font_size=22, color=GREEN,
            ),
        ).arrange(DOWN, buff=0.15, aligned_edge=LEFT)
        var_info.to_corner(DL).shift(UP * 0.5 + RIGHT * 0.3)
        var_box = SurroundingRectangle(var_info, color=WHITE, buff=0.12,
                                        fill_opacity=0.05)

        self.play(FadeIn(var_info), Create(var_box), run_time=1)
        self.wait(1)

        # ---------- 阶段7: 总结 ----------
        summary = Text(
            "Eigenvectors of covariance matrix = principal component directions",
            font_size=22, color=GREEN_A,
        ).to_edge(DOWN)
        self.play(FadeIn(summary), run_time=1)
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)
