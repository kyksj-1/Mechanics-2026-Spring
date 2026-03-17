"""
第3章：二次型 - Manim 可视化动画
===========================================

文件用途：为线性代数教学第3章制作3个可视化场景
输出格式：MP4 动画
运行命令：
  conda activate teaching_env && manim -qm ch3_quadratic_form_visualizations.py QuadraticFormSurface
  conda activate teaching_env && manim -qm ch3_quadratic_form_visualizations.py DiagonalizationOfQuadraticForm
  conda activate teaching_env && manim -qm ch3_quadratic_form_visualizations.py EnergyLandscape

作者：kyksj-1
"""

from manim import *
import numpy as np


# ============================================================
# Scene 1: 二次型曲面可视化
# 展示四种不同定性的二次型的三维曲面，配合等高线、矩阵和特征值标注
# ============================================================
class QuadraticFormSurface(ThreeDScene):
    def construct(self):
        # ---------- 定义四种二次型的参数 ----------
        cases = [
            {
                "name": "Positive Definite",
                "a": 2, "b": 0, "c": 3,
                "color": BLUE_C,
                "matrix_tex": r"A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}",
                "eigen_tex": r"\lambda_1 = 2 > 0,\; \lambda_2 = 3 > 0",
                "desc": "Bowl shape (elliptic paraboloid)",
                "z_range": [0, 6],
                "z_scale": 1.0,
            },
            {
                "name": "Negative Definite",
                "a": -2, "b": 0, "c": -3,
                "color": RED_C,
                "matrix_tex": r"A = \begin{pmatrix} -2 & 0 \\ 0 & -3 \end{pmatrix}",
                "eigen_tex": r"\lambda_1 = -2 < 0,\; \lambda_2 = -3 < 0",
                "desc": "Inverted bowl",
                "z_range": [-6, 0],
                "z_scale": 1.0,
            },
            {
                "name": "Indefinite",
                "a": 1, "b": 0, "c": -1,
                "color": ORANGE,
                "matrix_tex": r"A = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}",
                "eigen_tex": r"\lambda_1 = 1 > 0,\; \lambda_2 = -1 < 0",
                "desc": "Saddle surface (hyperbolic paraboloid)",
                "z_range": [-3, 3],
                "z_scale": 1.0,
            },
            {
                "name": "Positive Semi-Definite",
                "a": 1, "b": 0, "c": 0,
                "color": GREEN_C,
                "matrix_tex": r"A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}",
                "eigen_tex": r"\lambda_1 = 1 > 0,\; \lambda_2 = 0",
                "desc": "Cylindrical surface",
                "z_range": [0, 4],
                "z_scale": 1.0,
            },
        ]

        # ---------- 阶段1: 标题引入 ----------
        # 使用固定相机方式在平面上显示标题
        self.set_camera_orientation(phi=0, theta=-PI / 2)

        main_title = Text("Quadratic Form: Q(x,y) = ax^2 + 2bxy + cy^2",
                          font_size=30).to_edge(UP)
        formula = MathTex(
            r"Q(\vec{x}) = \vec{x}^T A \vec{x}",
            font_size=36,
        ).next_to(main_title, DOWN, buff=0.4)

        self.add_fixed_in_frame_mobjects(main_title, formula)
        self.play(Write(main_title), Write(formula), run_time=1.5)
        self.wait(1)
        self.play(FadeOut(main_title), FadeOut(formula), run_time=0.8)
        self.remove(main_title, formula)

        # ---------- 阶段2: 设置3D相机和坐标轴 ----------
        self.set_camera_orientation(phi=65 * DEGREES, theta=-45 * DEGREES)

        axes = ThreeDAxes(
            x_range=[-2.5, 2.5, 1],
            y_range=[-2.5, 2.5, 1],
            z_range=[-6, 6, 2],
            x_length=5,
            y_length=5,
            z_length=5,
            axis_config={"color": GREY_B, "stroke_width": 1.5},
        )
        # 坐标轴标签（使用 Text 避免 LaTeX 在 Windows 上的文件锁问题）
        x_label = axes.get_x_axis_label(Text("x", font_size=22))
        y_label = axes.get_y_axis_label(Text("y", font_size=22))
        z_label = axes.get_z_axis_label(Text("Q", font_size=22))

        self.play(Create(axes), Write(x_label), Write(y_label), Write(z_label),
                  run_time=1.5)
        self.wait(0.5)

        # ---------- 阶段3: 依次展示四种曲面 ----------
        prev_surface = None
        prev_contours = None
        prev_info = None

        for i, case in enumerate(cases):
            a, b, c = case["a"], case["b"], case["c"]

            # 创建三维曲面
            surface = Surface(
                lambda u, v, _a=a, _b=b, _c=c: axes.c2p(
                    u, v,
                    np.clip(_a * u**2 + 2 * _b * u * v + _c * v**2, -6, 6)
                ),
                u_range=[-2, 2],
                v_range=[-2, 2],
                resolution=(40, 40),
                fill_opacity=0.7,
                stroke_width=0.3,
                stroke_color=WHITE,
            )
            # 按z值设置颜色渐变
            surface.set_color_by_gradient(case["color"], WHITE, case["color"])

            # 等高线 (在z=z_min平面上投影)
            contour_group = VGroup()
            # 根据情形选择等高线数值
            if case["name"] == "Positive Definite":
                contour_levels = [0.5, 1.5, 3.0, 5.0]
            elif case["name"] == "Negative Definite":
                contour_levels = [-0.5, -1.5, -3.0, -5.0]
            elif case["name"] == "Indefinite":
                contour_levels = [-2, -0.5, 0.5, 2]
            else:
                contour_levels = [0.5, 1.5, 3.0]

            z_base = -6  # 等高线投影在底部平面
            for level in contour_levels:
                contour_pts = []
                num_pts = 200
                for t in np.linspace(0, 2 * np.pi, num_pts):
                    # 参数化等高线: a*x^2 + c*y^2 = level (b=0的情况简化)
                    if case["name"] == "Indefinite":
                        # 双曲线: x^2 - y^2 = level
                        if level > 0:
                            # 两支双曲线
                            for sign in [1, -1]:
                                pts = []
                                for s in np.linspace(-1.8, 1.8, 80):
                                    x_val = sign * np.sqrt(abs(level) + s**2)
                                    y_val = s
                                    if abs(x_val) <= 2:
                                        pts.append(axes.c2p(x_val, y_val, z_base))
                                if len(pts) > 2:
                                    curve = VMobject()
                                    curve.set_points_smoothly(pts)
                                    curve.set_stroke(case["color"], width=1.5, opacity=0.6)
                                    contour_group.add(curve)
                            # 跳过后面的通用逻辑
                            continue
                        elif level < 0:
                            for sign in [1, -1]:
                                pts = []
                                for s in np.linspace(-1.8, 1.8, 80):
                                    x_val = s
                                    y_val = sign * np.sqrt(abs(level) + s**2)
                                    if abs(y_val) <= 2:
                                        pts.append(axes.c2p(x_val, y_val, z_base))
                                if len(pts) > 2:
                                    curve = VMobject()
                                    curve.set_points_smoothly(pts)
                                    curve.set_stroke(case["color"], width=1.5, opacity=0.6)
                                    contour_group.add(curve)
                            continue
                    elif case["name"] == "Positive Semi-Definite":
                        # x^2 = level -> x = +/- sqrt(level)
                        if level > 0:
                            sq = np.sqrt(level)
                            for sign in [1, -1]:
                                pts = [axes.c2p(sign * sq, y_val, z_base)
                                       for y_val in np.linspace(-2, 2, 50)]
                                if len(pts) > 2:
                                    curve = VMobject()
                                    curve.set_points_smoothly(pts)
                                    curve.set_stroke(case["color"], width=1.5, opacity=0.6)
                                    contour_group.add(curve)
                        continue
                    else:
                        # 椭圆: a*x^2 + c*y^2 = level
                        if a * level <= 0 or c * level <= 0:
                            continue
                        rx = np.sqrt(abs(level / a))
                        ry = np.sqrt(abs(level / c))
                        x_val = rx * np.cos(t)
                        y_val = ry * np.sin(t)
                        if abs(x_val) <= 2 and abs(y_val) <= 2:
                            contour_pts.append(axes.c2p(x_val, y_val, z_base))

                if contour_pts and len(contour_pts) > 10:
                    contour_pts.append(contour_pts[0])  # 闭合曲线
                    curve = VMobject()
                    curve.set_points_smoothly(contour_pts)
                    curve.set_stroke(case["color"], width=1.5, opacity=0.6)
                    contour_group.add(curve)

            # 信息面板 (固定在画面前)
            info_title = Text(case["name"], font_size=28,
                              color=case["color"]).to_corner(UL).shift(DOWN * 0.3)
            info_matrix = MathTex(case["matrix_tex"],
                                  font_size=28).next_to(info_title, DOWN, buff=0.25)
            info_eigen = MathTex(case["eigen_tex"],
                                 font_size=24).next_to(info_matrix, DOWN, buff=0.2)
            info_desc = Text(case["desc"], font_size=20,
                             color=GREY_A).next_to(info_eigen, DOWN, buff=0.2)
            info_box = SurroundingRectangle(
                VGroup(info_title, info_matrix, info_eigen, info_desc),
                color=case["color"], buff=0.15, fill_opacity=0.08,
                stroke_width=1.5,
            )
            info_group = VGroup(info_title, info_matrix, info_eigen,
                                info_desc, info_box)

            # 将信息面板固定在画面前
            self.add_fixed_in_frame_mobjects(info_group)

            if i == 0:
                # 第一个曲面：正常创建
                self.play(Create(surface), run_time=2)
                self.play(Create(contour_group), run_time=1)
                self.play(FadeIn(info_group), run_time=1)
            else:
                # 后续曲面：变换切换
                if prev_info is not None:
                    self.play(FadeOut(prev_info), run_time=0.5)
                    self.remove(prev_info)
                self.play(
                    ReplacementTransform(prev_surface, surface),
                    ReplacementTransform(prev_contours, contour_group),
                    run_time=2,
                )
                self.play(FadeIn(info_group), run_time=0.8)

            # 慢慢旋转相机以展示3D效果
            self.begin_ambient_camera_rotation(rate=0.15)
            self.wait(3)
            self.stop_ambient_camera_rotation()

            prev_surface = surface
            prev_contours = contour_group
            prev_info = info_group

        # ---------- 阶段4: 总结 ----------
        self.play(FadeOut(prev_info), run_time=0.5)
        self.remove(prev_info)

        summary = Text(
            "Definiteness of A determines the shape of Q(x)",
            font_size=26, color=GREEN_A,
        ).to_edge(DOWN)
        self.add_fixed_in_frame_mobjects(summary)
        self.play(FadeIn(summary), run_time=1)
        self.wait(2)
        self.play(
            *[FadeOut(m) for m in self.mobjects],
            FadeOut(summary),
            run_time=1.5,
        )


# ============================================================
# Scene 2: 二次型的正交对角化
# 以 Q(x,y) = 5x^2 + 4xy + 2y^2 为例
# 矩阵 A = [[5,2],[2,2]]
# 特征值: lambda1 = 1, lambda2 = 6
# 特征向量: v1 = (-1,2)/sqrt5, v2 = (2,1)/sqrt5
# 旋转角 theta = arctan(2) ≈ 63.43 deg
# ============================================================
class DiagonalizationOfQuadraticForm(Scene):
    def construct(self):
        # ---------- 矩阵和特征分析 ----------
        A = np.array([[5, 2], [2, 2]])
        evals, evecs = np.linalg.eigh(A)
        # evals 按升序: lambda1=1, lambda2=6
        # evecs[:,0] 对应 lambda1=1, evecs[:,1] 对应 lambda2=6

        v1 = evecs[:, 0]  # 对应较小特征值
        v2 = evecs[:, 1]  # 对应较大特征值
        theta = np.arctan2(v1[1], v1[0])  # 旋转角

        # ---------- 阶段1: 标题和矩阵介绍 ----------
        title = Text("Diagonalization of Quadratic Form",
                      font_size=30).to_edge(UP)
        self.play(Write(title), run_time=1)

        # 展示二次型及其矩阵
        q_formula = MathTex(
            r"Q(x,y) = 5x^2 + 4xy + 2y^2",
            font_size=34,
        ).shift(UP * 2)
        q_matrix = MathTex(
            r"A = \begin{pmatrix} 5 & 2 \\ 2 & 2 \end{pmatrix}",
            font_size=34,
        ).next_to(q_formula, DOWN, buff=0.4)

        eigen_info = MathTex(
            r"\lambda_1 = 1,\; \lambda_2 = 6",
            font_size=30,
            color=YELLOW,
        ).next_to(q_matrix, DOWN, buff=0.4)

        self.play(Write(q_formula), run_time=1.2)
        self.play(Write(q_matrix), run_time=1)
        self.play(Write(eigen_info), run_time=1)
        self.wait(1)

        # 清理进入图形阶段
        self.play(
            FadeOut(title), FadeOut(q_formula),
            q_matrix.animate.scale(0.8).to_corner(UR).shift(DOWN * 0.3),
            eigen_info.animate.scale(0.85).to_corner(UR).shift(DOWN * 1.2),
            run_time=1,
        )

        # ---------- 阶段2: 原始等高线（倾斜椭圆） ----------
        plane = NumberPlane(
            x_range=[-4, 4, 1],
            y_range=[-3, 3, 1],
            x_length=8,
            y_length=6,
            background_line_style={
                "stroke_color": BLUE_E,
                "stroke_width": 0.8,
                "stroke_opacity": 0.3,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1.5},
        ).shift(DOWN * 0.2)
        x_lbl = plane.get_x_axis_label(MathTex("x", font_size=24))
        y_lbl = plane.get_y_axis_label(MathTex("y", font_size=24))

        self.play(Create(plane), Write(x_lbl), Write(y_lbl), run_time=1)

        # 绘制多条等高线
        contour_levels = [2, 5, 10, 18]
        contour_colors = [BLUE_A, BLUE_C, BLUE_D, BLUE_E]
        original_contours = VGroup()

        for level, col in zip(contour_levels, contour_colors):
            contour_pts = []
            for t in np.linspace(0, 2 * np.pi, 200):
                # 参数化: (x,y) = R(theta) @ (sqrt(level/lam1)*cos(t), sqrt(level/lam2)*sin(t))
                u = np.sqrt(level / evals[0]) * np.cos(t)
                w = np.sqrt(level / evals[1]) * np.sin(t)
                # 旋转回原坐标: [x,y] = P @ [u,w]
                xy = evecs @ np.array([u, w])
                if abs(xy[0]) <= 3.8 and abs(xy[1]) <= 2.8:
                    contour_pts.append(plane.c2p(xy[0], xy[1]))
            if len(contour_pts) > 10:
                contour_pts.append(contour_pts[0])
                curve = VMobject()
                curve.set_points_smoothly(contour_pts)
                curve.set_stroke(col, width=2, opacity=0.8)
                original_contours.add(curve)

        self.play(
            LaggedStart(*[Create(c) for c in original_contours], lag_ratio=0.2),
            run_time=2,
        )

        original_label = Text("Original ellipses (tilted)",
                              font_size=22, color=GREY_A).to_edge(DOWN)
        self.play(FadeIn(original_label), run_time=0.5)
        self.wait(1.5)

        # ---------- 阶段3: 标注特征向量方向 ----------
        # 特征向量1 (对应 lambda1)
        ev_arrow1 = Arrow(
            plane.c2p(0, 0), plane.c2p(*(v1 * 2.5)),
            buff=0, stroke_width=4, color=YELLOW,
            max_tip_length_to_length_ratio=0.12,
        )
        ev_arrow2 = Arrow(
            plane.c2p(0, 0), plane.c2p(*(v2 * 2.5)),
            buff=0, stroke_width=4, color=GREEN_B,
            max_tip_length_to_length_ratio=0.12,
        )
        ev1_lbl = MathTex(r"\vec{v}_1\;(\lambda_1=1)", font_size=22,
                          color=YELLOW).next_to(ev_arrow1.get_end(), UR, buff=0.1)
        ev2_lbl = MathTex(r"\vec{v}_2\;(\lambda_2=6)", font_size=22,
                          color=GREEN_B).next_to(ev_arrow2.get_end(), UL, buff=0.1)

        self.play(FadeOut(original_label), run_time=0.3)
        self.play(
            GrowArrow(ev_arrow1), GrowArrow(ev_arrow2),
            Write(ev1_lbl), Write(ev2_lbl),
            run_time=1.5,
        )

        # 标注旋转角
        angle_arc = Arc(
            radius=0.8,
            start_angle=0,
            angle=theta,
            color=TEAL,
            stroke_width=2,
        ).shift(plane.c2p(0, 0))
        angle_label = MathTex(r"\theta", font_size=22,
                              color=TEAL).next_to(angle_arc, RIGHT, buff=0.1)

        self.play(Create(angle_arc), Write(angle_label), run_time=0.8)
        self.wait(1.5)

        # ---------- 阶段4: 旋转坐标轴对齐特征方向 ----------
        rotate_text = Text("Rotate axes to align with eigenvectors",
                           font_size=22, color=TEAL_A).to_edge(DOWN)
        self.play(FadeIn(rotate_text), run_time=0.5)

        # 创建新坐标轴（沿特征方向）
        new_x_axis = Arrow(
            plane.c2p(*(- v1 * 3.5)), plane.c2p(*(v1 * 3.5)),
            buff=0, stroke_width=2.5, color=YELLOW,
            max_tip_length_to_length_ratio=0.06,
        )
        new_y_axis = Arrow(
            plane.c2p(*(- v2 * 2.8)), plane.c2p(*(v2 * 2.8)),
            buff=0, stroke_width=2.5, color=GREEN_B,
            max_tip_length_to_length_ratio=0.06,
        )
        new_x_lbl = MathTex(r"x'", font_size=24,
                            color=YELLOW).next_to(new_x_axis.get_end(), DR, buff=0.1)
        new_y_lbl = MathTex(r"y'", font_size=24,
                            color=GREEN_B).next_to(new_y_axis.get_end(), UL, buff=0.1)

        self.play(
            Create(new_x_axis), Create(new_y_axis),
            Write(new_x_lbl), Write(new_y_lbl),
            run_time=1.5,
        )
        self.wait(1)

        # ---------- 阶段5: 变换到标准形（旋转一切使特征方向与坐标轴对齐） ----------
        self.play(FadeOut(rotate_text), run_time=0.3)

        # 计算旋转矩阵
        cos_t = np.cos(-theta)
        sin_t = np.sin(-theta)

        # 创建变换后的等高线（轴对齐的椭圆）
        standard_contours = VGroup()
        for level, col in zip(contour_levels, contour_colors):
            contour_pts = []
            rx = np.sqrt(level / evals[0])
            ry = np.sqrt(level / evals[1])
            for t in np.linspace(0, 2 * np.pi, 200):
                x_val = rx * np.cos(t)
                y_val = ry * np.sin(t)
                if abs(x_val) <= 3.8 and abs(y_val) <= 2.8:
                    contour_pts.append(plane.c2p(x_val, y_val))
            if len(contour_pts) > 10:
                contour_pts.append(contour_pts[0])
                curve = VMobject()
                curve.set_points_smoothly(contour_pts)
                curve.set_stroke(col, width=2, opacity=0.8)
                standard_contours.add(curve)

        # 新的标准坐标轴标签
        new_plane_x_lbl = MathTex(r"x'", font_size=24,
                                  color=YELLOW).next_to(plane.get_x_axis().get_end(), DR, buff=0.1)
        new_plane_y_lbl = MathTex(r"y'", font_size=24,
                                  color=GREEN_B).next_to(plane.get_y_axis().get_end(), UL, buff=0.1)

        transform_text = Text("Standard form: axes aligned with eigenvectors",
                              font_size=22, color=GREEN_A).to_edge(DOWN)

        # 执行旋转变换
        self.play(
            ReplacementTransform(original_contours, standard_contours),
            FadeOut(ev_arrow1), FadeOut(ev_arrow2),
            FadeOut(ev1_lbl), FadeOut(ev2_lbl),
            FadeOut(angle_arc), FadeOut(angle_label),
            FadeOut(new_x_axis), FadeOut(new_y_axis),
            FadeOut(new_x_lbl), FadeOut(new_y_lbl),
            FadeOut(x_lbl), FadeOut(y_lbl),
            FadeIn(new_plane_x_lbl), FadeIn(new_plane_y_lbl),
            FadeIn(transform_text),
            run_time=2.5,
        )
        self.wait(1)

        # ---------- 阶段6: 标注标准形公式 ----------
        standard_form = MathTex(
            r"Q = \lambda_1 {x'}^2 + \lambda_2 {y'}^2 = {x'}^2 + 6{y'}^2",
            font_size=32,
            color=YELLOW,
        ).to_corner(UL).shift(DOWN * 0.5 + RIGHT * 0.2)
        std_box = SurroundingRectangle(standard_form, color=YELLOW, buff=0.12,
                                       fill_opacity=0.08)

        relation = MathTex(
            r"\begin{pmatrix} x \\ y \end{pmatrix} = P "
            r"\begin{pmatrix} x' \\ y' \end{pmatrix}",
            font_size=28,
        ).next_to(standard_form, DOWN, buff=0.3)

        self.play(
            Write(standard_form), Create(std_box),
            run_time=1.2,
        )
        self.play(Write(relation), run_time=1)
        self.wait(1.5)

        # ---------- 阶段7: 总结文字 ----------
        self.play(FadeOut(transform_text), run_time=0.3)
        summary = Text(
            "Orthogonal diagonalization aligns axes with eigenvectors",
            font_size=24, color=GREEN_A,
        ).to_edge(DOWN)
        self.play(FadeIn(summary), run_time=1)
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 3: 能量地形
# 展示二次型作为能量面的概念：等高线、梯度场、粒子轨迹
# ============================================================
class EnergyLandscape(Scene):
    def construct(self):
        # ---------- 阶段1: 标题 ----------
        title = Text("Quadratic Form as Energy Landscape",
                      font_size=30).to_edge(UP)
        self.play(Write(title), run_time=1)
        self.wait(0.5)

        # ===== Part A: 正定情形 - 稳定平衡 =====
        self._show_positive_definite(title)

        # ===== Part B: 不定情形 - 鞍点 =====
        self._show_indefinite(title)

        self.wait(1)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)

    def _show_positive_definite(self, title):
        """正定情形：梯度指向唯一极小值，粒子滚向最低点"""
        A = np.array([[2, 0.5], [0.5, 1.5]])

        # 标题
        case_title = Text("Positive Definite: Stable Minimum",
                          font_size=24, color=BLUE_C).next_to(title, DOWN, buff=0.3)
        mat_label = MathTex(
            r"A = \begin{pmatrix} 2 & 0.5 \\ 0.5 & 1.5 \end{pmatrix}",
            font_size=26,
        ).to_corner(UR).shift(DOWN * 0.6 + LEFT * 0.2)
        mat_bg = SurroundingRectangle(mat_label, color=BLUE_C, buff=0.1,
                                       fill_opacity=0.05)

        self.play(Write(case_title), Write(mat_label), Create(mat_bg), run_time=1)

        # 坐标系
        plane = NumberPlane(
            x_range=[-3.5, 3.5, 1],
            y_range=[-2.5, 2.5, 1],
            x_length=7,
            y_length=5,
            background_line_style={
                "stroke_color": BLUE_E,
                "stroke_width": 0.6,
                "stroke_opacity": 0.2,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1.2},
        ).shift(DOWN * 0.3)

        self.play(Create(plane), run_time=0.8)

        # 等高线
        evals_pd, evecs_pd = np.linalg.eigh(A)
        contour_levels = [0.5, 1.5, 3, 5, 8]
        contour_colors = color_gradient([BLUE_A, BLUE_D], len(contour_levels))
        contours = VGroup()

        for level, col in zip(contour_levels, contour_colors):
            pts = []
            for t in np.linspace(0, 2 * np.pi, 200):
                u = np.sqrt(level / evals_pd[0]) * np.cos(t)
                w = np.sqrt(level / evals_pd[1]) * np.sin(t)
                xy = evecs_pd @ np.array([u, w])
                if abs(xy[0]) <= 3.3 and abs(xy[1]) <= 2.3:
                    pts.append(plane.c2p(xy[0], xy[1]))
            if len(pts) > 10:
                pts.append(pts[0])
                curve = VMobject()
                curve.set_points_smoothly(pts)
                curve.set_stroke(col, width=1.5, opacity=0.7)
                contours.add(curve)

        self.play(
            LaggedStart(*[Create(c) for c in contours], lag_ratio=0.15),
            run_time=1.5,
        )

        # 梯度场 (负梯度方向，即下降方向)
        # 梯度 grad Q = 2Ax
        def neg_grad_func(pos):
            coords = np.array(plane.p2c(pos))
            grad = 2 * A @ coords
            # 负梯度（下降方向）
            neg_grad = -grad
            speed = np.linalg.norm(neg_grad)
            if speed < 0.01:
                return np.array([0, 0, 0])
            scale = min(0.35, 0.35 * speed / 2.0)
            return np.array([
                neg_grad[0] * scale * plane.get_x_unit_size(),
                neg_grad[1] * scale * plane.get_y_unit_size(),
                0
            ])

        # 梯度箭头场
        arrow_field = ArrowVectorField(
            neg_grad_func,
            x_range=[plane.c2p(-3, 0)[0], plane.c2p(3, 0)[0], 0.7],
            y_range=[plane.c2p(0, -2)[1], plane.c2p(0, 2)[1], 0.7],
            length_func=lambda norm: min(norm, 0.4),
            color=GREY_A,
            opacity=0.5,
        )

        self.play(Create(arrow_field), run_time=1.5)

        # 极小值点标注
        min_dot = Dot(plane.c2p(0, 0), color=RED, radius=0.1)
        min_label = Text("Minimum", font_size=18,
                         color=RED).next_to(min_dot, DOWN, buff=0.15)
        self.play(FadeIn(min_dot, scale=1.5), Write(min_label), run_time=0.8)
        self.wait(0.5)

        # 粒子轨迹动画：从不同初始位置"滚下"
        start_points = [
            np.array([2.5, 1.5]),
            np.array([-2.0, -1.5]),
            np.array([1.0, -2.0]),
            np.array([-2.5, 1.0]),
        ]
        traj_colors = [YELLOW, GREEN_B, TEAL, ORANGE]

        desc_text = Text("Particles roll toward the minimum",
                         font_size=20, color=GREY_A).to_edge(DOWN)
        self.play(FadeIn(desc_text), run_time=0.5)

        # 预计算轨迹
        dt = 0.02
        num_steps = 200
        trajectories = []
        for sp in start_points:
            traj = [sp.copy()]
            pos = sp.copy()
            for _ in range(num_steps):
                grad = 2 * A @ pos
                pos = pos - dt * grad
                traj.append(pos.copy())
            trajectories.append(traj)

        # 动画：同时展示多条轨迹
        traj_lines = VGroup()
        particles = VGroup()

        for traj, col in zip(trajectories, traj_colors):
            # 绘制轨迹线
            traj_pts = [plane.c2p(p[0], p[1]) for p in traj
                        if abs(p[0]) <= 3.5 and abs(p[1]) <= 2.5]
            if len(traj_pts) > 2:
                line = VMobject()
                line.set_points_smoothly(traj_pts)
                line.set_stroke(col, width=2.5, opacity=0.8)
                traj_lines.add(line)

            # 起始粒子
            dot = Dot(plane.c2p(traj[0][0], traj[0][1]),
                      color=col, radius=0.08)
            particles.add(dot)

        # 先显示起始点
        self.play(
            LaggedStart(*[FadeIn(p, scale=1.5) for p in particles], lag_ratio=0.1),
            run_time=0.8,
        )

        # 展示轨迹
        self.play(
            LaggedStart(*[Create(l) for l in traj_lines], lag_ratio=0.1),
            run_time=3,
        )

        # 粒子移动到终点
        for i, (particle, traj) in enumerate(zip(particles, trajectories)):
            end_pos = traj[-1]
            particle.generate_target()
            particle.target.move_to(plane.c2p(end_pos[0], end_pos[1]))

        self.play(
            *[MoveToTarget(p) for p in particles],
            run_time=1,
        )
        self.wait(2)

        # 清除正定情形
        self.play(
            FadeOut(case_title), FadeOut(mat_label), FadeOut(mat_bg),
            FadeOut(plane), FadeOut(contours), FadeOut(arrow_field),
            FadeOut(min_dot), FadeOut(min_label), FadeOut(desc_text),
            FadeOut(traj_lines), FadeOut(particles),
            run_time=1,
        )

    def _show_indefinite(self, title):
        """不定情形：鞍点，梯度为零但非极值"""
        A = np.array([[1, 0], [0, -1]])

        case_title = Text("Indefinite: Saddle Point",
                          font_size=24, color=ORANGE).next_to(title, DOWN, buff=0.3)
        mat_label = MathTex(
            r"A = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}",
            font_size=26,
        ).to_corner(UR).shift(DOWN * 0.6 + LEFT * 0.2)
        mat_bg = SurroundingRectangle(mat_label, color=ORANGE, buff=0.1,
                                       fill_opacity=0.05)

        self.play(Write(case_title), Write(mat_label), Create(mat_bg), run_time=1)

        # 坐标系
        plane = NumberPlane(
            x_range=[-3.5, 3.5, 1],
            y_range=[-2.5, 2.5, 1],
            x_length=7,
            y_length=5,
            background_line_style={
                "stroke_color": ORANGE,
                "stroke_width": 0.5,
                "stroke_opacity": 0.15,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1.2},
        ).shift(DOWN * 0.3)

        self.play(Create(plane), run_time=0.8)

        # 等高线: x^2 - y^2 = c (双曲线)
        contour_levels = [-3, -1.5, -0.5, 0.5, 1.5, 3]
        contours = VGroup()

        for level in contour_levels:
            if abs(level) < 0.01:
                continue
            if level > 0:
                # x^2 - y^2 = level > 0: x = +/- sqrt(level + y^2)
                color = RED_C
                for sign in [1, -1]:
                    pts = []
                    for s in np.linspace(-2.2, 2.2, 120):
                        x_val = sign * np.sqrt(level + s**2)
                        y_val = s
                        if abs(x_val) <= 3.3 and abs(y_val) <= 2.3:
                            pts.append(plane.c2p(x_val, y_val))
                    if len(pts) > 5:
                        curve = VMobject()
                        curve.set_points_smoothly(pts)
                        curve.set_stroke(color, width=1.5, opacity=0.7)
                        contours.add(curve)
            else:
                # x^2 - y^2 = level < 0: y = +/- sqrt(-level + x^2)
                color = BLUE_C
                for sign in [1, -1]:
                    pts = []
                    for s in np.linspace(-3.0, 3.0, 120):
                        x_val = s
                        y_sq = -level + s**2
                        if y_sq >= 0:
                            y_val = sign * np.sqrt(y_sq)
                            if abs(x_val) <= 3.3 and abs(y_val) <= 2.3:
                                pts.append(plane.c2p(x_val, y_val))
                    if len(pts) > 5:
                        curve = VMobject()
                        curve.set_points_smoothly(pts)
                        curve.set_stroke(color, width=1.5, opacity=0.7)
                        contours.add(curve)

        # 渐近线 (y = x 和 y = -x)
        asymp1 = DashedLine(
            plane.c2p(-2.5, -2.5), plane.c2p(2.5, 2.5),
            color=GREY_A, stroke_width=1.5,
        )
        asymp2 = DashedLine(
            plane.c2p(-2.5, 2.5), plane.c2p(2.5, -2.5),
            color=GREY_A, stroke_width=1.5,
        )

        self.play(
            Create(asymp1), Create(asymp2),
            LaggedStart(*[Create(c) for c in contours], lag_ratio=0.1),
            run_time=2,
        )

        # 梯度场
        def neg_grad_func_indef(pos):
            coords = np.array(plane.p2c(pos))
            grad = 2 * A @ coords
            neg_grad = -grad
            speed = np.linalg.norm(neg_grad)
            if speed < 0.01:
                return np.array([0, 0, 0])
            scale = min(0.35, 0.35 * speed / 2.0)
            return np.array([
                neg_grad[0] * scale * plane.get_x_unit_size(),
                neg_grad[1] * scale * plane.get_y_unit_size(),
                0
            ])

        arrow_field = ArrowVectorField(
            neg_grad_func_indef,
            x_range=[plane.c2p(-3, 0)[0], plane.c2p(3, 0)[0], 0.7],
            y_range=[plane.c2p(0, -2)[1], plane.c2p(0, 2)[1], 0.7],
            length_func=lambda norm: min(norm, 0.4),
            color=GREY_A,
            opacity=0.4,
        )
        self.play(Create(arrow_field), run_time=1.5)

        # 鞍点标注
        saddle_dot = Dot(plane.c2p(0, 0), color=ORANGE, radius=0.12)
        saddle_label = Text("Saddle Point", font_size=18,
                            color=ORANGE).next_to(saddle_dot, DL, buff=0.15)
        grad_zero = MathTex(r"\nabla Q = 0", font_size=22,
                            color=ORANGE).next_to(saddle_label, DOWN, buff=0.1)

        self.play(
            FadeIn(saddle_dot, scale=1.5),
            Write(saddle_label), Write(grad_zero),
            run_time=1,
        )
        self.wait(1)

        # 粒子轨迹：有些逃逸，有些趋近再逃逸
        start_points = [
            np.array([0.3, 2.0]),     # 从上方开始, y方向不稳定
            np.array([-0.3, -2.0]),    # 从下方开始
            np.array([2.0, 0.3]),     # 从右方开始, x方向稳定
            np.array([-2.0, -0.3]),    # 从左方开始
        ]
        traj_colors = [GREEN_B, GREEN_B, YELLOW, YELLOW]

        desc_text = Text("Some directions attract, others repel",
                         font_size=20, color=GREY_A).to_edge(DOWN)
        self.play(FadeIn(desc_text), run_time=0.5)

        # 预计算轨迹
        dt = 0.015
        num_steps = 250
        trajectories = []
        for sp in start_points:
            traj = [sp.copy()]
            pos = sp.copy()
            for _ in range(num_steps):
                grad = 2 * A @ pos
                pos = pos - dt * grad
                # 限制范围
                pos = np.clip(pos, -3.5, 3.5)
                traj.append(pos.copy())
            trajectories.append(traj)

        # 绘制轨迹
        traj_lines = VGroup()
        particles = VGroup()

        for traj, col in zip(trajectories, traj_colors):
            traj_pts = [plane.c2p(p[0], p[1]) for p in traj
                        if abs(p[0]) <= 3.5 and abs(p[1]) <= 2.5]
            if len(traj_pts) > 2:
                line = VMobject()
                line.set_points_smoothly(traj_pts)
                line.set_stroke(col, width=2.5, opacity=0.8)
                traj_lines.add(line)

            dot = Dot(plane.c2p(traj[0][0], traj[0][1]),
                      color=col, radius=0.08)
            particles.add(dot)

        self.play(
            LaggedStart(*[FadeIn(p, scale=1.5) for p in particles], lag_ratio=0.1),
            run_time=0.6,
        )
        self.play(
            LaggedStart(*[Create(l) for l in traj_lines], lag_ratio=0.1),
            run_time=3,
        )

        # 粒子移到终点
        for particle, traj in zip(particles, trajectories):
            end_pos = traj[-1]
            particle.generate_target()
            particle.target.move_to(plane.c2p(
                np.clip(end_pos[0], -3.3, 3.3),
                np.clip(end_pos[1], -2.3, 2.3)
            ))

        self.play(*[MoveToTarget(p) for p in particles], run_time=1)
        self.wait(1)

        # 说明标签
        x_dir_label = Text("x-direction: stable (descends)",
                           font_size=18, color=YELLOW).to_corner(DL).shift(UP * 0.8)
        y_dir_label = Text("y-direction: unstable (ascends)",
                           font_size=18, color=GREEN_B).next_to(x_dir_label, DOWN, buff=0.15)
        self.play(FadeIn(x_dir_label), FadeIn(y_dir_label), run_time=0.8)
        self.wait(2)

        # 总结
        summary = Text(
            "Definiteness determines stability of critical points",
            font_size=24, color=GREEN_A,
        ).to_edge(DOWN)
        self.play(
            FadeOut(desc_text),
            FadeIn(summary),
            run_time=1,
        )
        self.wait(2)
