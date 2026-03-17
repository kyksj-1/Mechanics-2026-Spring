"""
第1章：特征值与特征向量 - Manim 可视化动画
===========================================

文件用途：为线性代数教学第1章制作4个可视化场景
输出格式：MP4 动画
运行命令：
  conda activate teaching_env && manim -qm ch1_eigenvalue_visualizations.py EigenVectorGeometry
  conda activate teaching_env && manim -qm ch1_eigenvalue_visualizations.py CharacteristicPolynomial
  conda activate teaching_env && manim -qm ch1_eigenvalue_visualizations.py EigenvalueSpectrum
  conda activate teaching_env && manim -qm ch1_eigenvalue_visualizations.py LinearDynamics

作者：kyksj-1
"""

from manim import *
import numpy as np


# ============================================================
# Scene 1: 特征向量的几何直觉
# 展示矩阵 A = [[2,1],[1,2]] 对平面向量的变换效果
# 特征值: lambda1=3 (方向(1,1)), lambda2=1 (方向(1,-1))
# ============================================================
class EigenVectorGeometry(Scene):
    def construct(self):
        # ---------- 阶段1: 搭建坐标平面和初始向量 ----------
        plane = NumberPlane(
            x_range=[-5, 5, 1],
            y_range=[-4, 4, 1],
            background_line_style={
                "stroke_color": BLUE_E,
                "stroke_width": 1,
                "stroke_opacity": 0.4,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 2},
        )
        plane.scale(0.85)

        # 矩阵 A 的定义
        A = np.array([[2, 1], [1, 2]])

        # 标题
        title = Text("Linear Transformation by Matrix A", font_size=32).to_edge(UP)
        matrix_label = MathTex(
            r"A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}",
            font_size=36,
        ).to_corner(UR).shift(DOWN * 0.5 + LEFT * 0.3)
        matrix_bg = SurroundingRectangle(
            matrix_label, color=BLUE, buff=0.15, fill_opacity=0.1
        )

        self.play(Create(plane), run_time=1.5)
        self.play(
            Write(title),
            FadeIn(matrix_label),
            Create(matrix_bg),
            run_time=1.2,
        )
        self.wait(0.5)

        # ---------- 阶段2: 绘制均匀分布的向量 ----------
        num_vectors = 12
        angles = np.linspace(0, 2 * np.pi, num_vectors, endpoint=False)
        radius = 2.0
        # 渐变色从蓝到绿
        colors = color_gradient([BLUE_C, GREEN_C, TEAL_C], num_vectors)

        original_vectors = VGroup()
        for i, angle in enumerate(angles):
            end_point = plane.c2p(
                radius * np.cos(angle), radius * np.sin(angle)
            )
            vec = Arrow(
                plane.c2p(0, 0),
                end_point,
                buff=0,
                stroke_width=3,
                color=colors[i],
                max_tip_length_to_length_ratio=0.15,
            )
            original_vectors.add(vec)

        before_label = Text(
            "Before transformation", font_size=24, color=GREY_A
        ).to_edge(DOWN)

        self.play(
            LaggedStart(
                *[GrowArrow(v) for v in original_vectors],
                lag_ratio=0.08,
            ),
            FadeIn(before_label),
            run_time=2,
        )
        self.wait(1)

        # ---------- 阶段3: 施加线性变换 ----------
        # 计算变换后的向量终点
        transformed_vectors = VGroup()
        for i, angle in enumerate(angles):
            original_tip = np.array([radius * np.cos(angle), radius * np.sin(angle)])
            new_tip = A @ original_tip
            end_point = plane.c2p(new_tip[0], new_tip[1])
            vec = Arrow(
                plane.c2p(0, 0),
                end_point,
                buff=0,
                stroke_width=3,
                color=colors[i],
                max_tip_length_to_length_ratio=0.15,
            )
            transformed_vectors.add(vec)

        after_label = Text(
            "After transformation: v -> Av", font_size=24, color=GREY_A
        ).to_edge(DOWN)

        self.play(
            FadeOut(before_label),
            run_time=0.3,
        )
        self.play(
            *[
                Transform(original_vectors[i], transformed_vectors[i])
                for i in range(num_vectors)
            ],
            FadeIn(after_label),
            run_time=2.5,
        )
        self.wait(1)

        # ---------- 阶段4: 高亮特征向量方向 ----------
        # 特征向量1: (1,1), 特征值3
        ev1_dir = np.array([1, 1]) / np.sqrt(2)
        ev1_original_end = plane.c2p(*(ev1_dir * radius))
        ev1_transformed_end = plane.c2p(*(ev1_dir * radius * 3))

        eigen_arrow_1_before = Arrow(
            plane.c2p(0, 0),
            ev1_original_end,
            buff=0,
            stroke_width=5,
            color=YELLOW,
            max_tip_length_to_length_ratio=0.12,
        )
        eigen_arrow_1_after = Arrow(
            plane.c2p(0, 0),
            ev1_transformed_end,
            buff=0,
            stroke_width=5,
            color=YELLOW,
            max_tip_length_to_length_ratio=0.12,
        )

        # 特征向量2: (1,-1), 特征值1
        ev2_dir = np.array([1, -1]) / np.sqrt(2)
        ev2_original_end = plane.c2p(*(ev2_dir * radius))
        ev2_transformed_end = plane.c2p(*(ev2_dir * radius * 1))

        eigen_arrow_2_before = Arrow(
            plane.c2p(0, 0),
            ev2_original_end,
            buff=0,
            stroke_width=5,
            color=RED_B,
            max_tip_length_to_length_ratio=0.15,
        )
        eigen_arrow_2_after = Arrow(
            plane.c2p(0, 0),
            ev2_transformed_end,
            buff=0,
            stroke_width=5,
            color=RED_B,
            max_tip_length_to_length_ratio=0.15,
        )

        # 先显示变换前的特征向量，再变换
        self.play(
            FadeOut(after_label),
            *[FadeOut(original_vectors[i]) for i in range(num_vectors)],
            run_time=0.8,
        )

        # 重新绘制所有向量（变换前状态）
        original_vectors_2 = VGroup()
        for i, angle in enumerate(angles):
            end_point = plane.c2p(
                radius * np.cos(angle), radius * np.sin(angle)
            )
            vec = Arrow(
                plane.c2p(0, 0),
                end_point,
                buff=0,
                stroke_width=2.5,
                color=colors[i],
                max_tip_length_to_length_ratio=0.15,
            ).set_opacity(0.4)
            original_vectors_2.add(vec)

        self.play(
            *[FadeIn(v) for v in original_vectors_2],
            GrowArrow(eigen_arrow_1_before),
            GrowArrow(eigen_arrow_2_before),
            run_time=1.2,
        )

        # 标注特征向量方向
        ev1_label = MathTex(
            r"\vec{v}_1 = \begin{pmatrix}1\\1\end{pmatrix}",
            font_size=28,
            color=YELLOW,
        ).next_to(eigen_arrow_1_before.get_end(), UR, buff=0.15)

        ev2_label = MathTex(
            r"\vec{v}_2 = \begin{pmatrix}1\\-1\end{pmatrix}",
            font_size=28,
            color=RED_B,
        ).next_to(eigen_arrow_2_before.get_end(), DR, buff=0.15)

        self.play(Write(ev1_label), Write(ev2_label), run_time=1)
        self.wait(1)

        # 施加变换：普通向量方向改变，特征向量只被拉伸
        transformed_vectors_2 = VGroup()
        for i, angle in enumerate(angles):
            original_tip = np.array([radius * np.cos(angle), radius * np.sin(angle)])
            new_tip = A @ original_tip
            end_point = plane.c2p(new_tip[0], new_tip[1])
            vec = Arrow(
                plane.c2p(0, 0),
                end_point,
                buff=0,
                stroke_width=2.5,
                color=colors[i],
                max_tip_length_to_length_ratio=0.15,
            ).set_opacity(0.4)
            transformed_vectors_2.add(vec)

        ev1_label_new = MathTex(
            r"A\vec{v}_1 = 3\vec{v}_1",
            font_size=28,
            color=YELLOW,
        ).next_to(eigen_arrow_1_after.get_end(), UR, buff=0.15)

        ev2_label_new = MathTex(
            r"A\vec{v}_2 = 1 \cdot \vec{v}_2",
            font_size=28,
            color=RED_B,
        ).next_to(eigen_arrow_2_after.get_end(), DR, buff=0.15)

        apply_label = Text(
            "Apply A", font_size=28, color=WHITE
        ).to_edge(DOWN)

        self.play(FadeIn(apply_label), run_time=0.5)
        self.play(
            *[
                Transform(original_vectors_2[i], transformed_vectors_2[i])
                for i in range(num_vectors)
            ],
            Transform(eigen_arrow_1_before, eigen_arrow_1_after),
            Transform(eigen_arrow_2_before, eigen_arrow_2_after),
            Transform(ev1_label, ev1_label_new),
            Transform(ev2_label, ev2_label_new),
            run_time=2.5,
        )
        self.wait(1)

        # ---------- 阶段5: 总结文字 ----------
        self.play(FadeOut(apply_label), run_time=0.3)
        summary = VGroup(
            Text("Eigenvector: direction preserved under transformation",
                 font_size=26, color=GREEN_A),
            Text("Other vectors: both direction and magnitude change",
                 font_size=26, color=GREY_B),
        ).arrange(DOWN, buff=0.3).to_edge(DOWN)

        self.play(FadeIn(summary), run_time=1)
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 2: 特征多项式求解过程可视化
# 从矩阵 A 出发，推导 det(A - lambda I) = 0
# 画出特征多项式图像，标注零点
# ============================================================
class CharacteristicPolynomial(Scene):
    def construct(self):
        # ---------- 阶段1: 展示矩阵 A ----------
        title = Text("Finding Eigenvalues: Characteristic Polynomial",
                      font_size=30).to_edge(UP)
        self.play(Write(title), run_time=1)

        mat_A = MathTex(
            r"A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}",
            font_size=42,
        )
        mat_A.shift(UP * 1.5)
        self.play(Write(mat_A), run_time=1.2)
        self.wait(0.8)

        # ---------- 阶段2: 构建 A - lambda I ----------
        step1_text = Text("Step 1: Form (A - lambda I)", font_size=24,
                          color=GREY_A).next_to(title, DOWN, buff=0.3)
        # 用 lambda 写法避免 LaTeX 冲突
        mat_AminusLI = MathTex(
            r"A - \lambda I = \begin{pmatrix} 2 - \lambda & 1 \\ 1 & 2 - \lambda \end{pmatrix}",
            font_size=40,
        )
        mat_AminusLI.shift(UP * 1.5)

        self.play(Write(step1_text), run_time=0.8)
        self.play(
            TransformMatchingTex(mat_A, mat_AminusLI),
            run_time=2,
        )
        self.wait(1)

        # ---------- 阶段3: 计算行列式 ----------
        step2_text = Text("Step 2: Compute det(A - lambda I) = 0",
                          font_size=24, color=GREY_A)
        step2_text.next_to(step1_text, DOWN, buff=0.2)

        det_expand = MathTex(
            r"\det(A - \lambda I) = (2-\lambda)^2 - 1",
            font_size=38,
        ).next_to(mat_AminusLI, DOWN, buff=0.6)

        det_simplify = MathTex(
            r"= \lambda^2 - 4\lambda + 3",
            font_size=38,
        ).next_to(det_expand, DOWN, buff=0.3)

        det_factor = MathTex(
            r"= (\lambda - 1)(\lambda - 3) = 0",
            font_size=38,
            color=YELLOW,
        ).next_to(det_simplify, DOWN, buff=0.3)

        self.play(Write(step2_text), run_time=0.8)
        self.play(Write(det_expand), run_time=1.5)
        self.wait(0.5)
        self.play(Write(det_simplify), run_time=1.2)
        self.wait(0.5)
        self.play(Write(det_factor), run_time=1.2)
        self.wait(1)

        # ---------- 阶段4: 清理并画函数图像 ----------
        eigenvalues_text = MathTex(
            r"\lambda_1 = 1, \quad \lambda_2 = 3",
            font_size=36,
            color=YELLOW,
        ).to_corner(UR).shift(DOWN * 0.8 + LEFT * 0.3)
        ev_box = SurroundingRectangle(eigenvalues_text, color=YELLOW, buff=0.15)

        self.play(
            FadeOut(title), FadeOut(step1_text), FadeOut(step2_text),
            FadeOut(mat_AminusLI), FadeOut(det_expand), FadeOut(det_simplify),
            FadeOut(det_factor),
            run_time=1,
        )

        # 函数图像
        axes = Axes(
            x_range=[-0.5, 5, 1],
            y_range=[-1.5, 4, 1],
            x_length=8,
            y_length=5,
            axis_config={"color": GREY_B, "include_numbers": True},
            tips=True,
        ).shift(DOWN * 0.3)
        x_label = axes.get_x_axis_label(MathTex(r"\lambda", font_size=32))
        y_label = axes.get_y_axis_label(
            MathTex(r"f(\lambda)", font_size=32), edge=LEFT, direction=LEFT
        )

        # 特征多项式 f(lambda) = lambda^2 - 4*lambda + 3
        graph = axes.plot(
            lambda x: x**2 - 4 * x + 3,
            x_range=[-0.3, 4.5],
            color=BLUE_C,
            stroke_width=3,
        )
        graph_label = MathTex(
            r"f(\lambda) = \lambda^2 - 4\lambda + 3",
            font_size=30,
            color=BLUE_C,
        ).next_to(graph, UP, buff=0.1).shift(RIGHT * 1)

        chart_title = Text(
            "Characteristic Polynomial", font_size=28
        ).to_edge(UP)

        self.play(
            Write(chart_title),
            Create(axes),
            Write(x_label),
            Write(y_label),
            run_time=1.5,
        )
        self.play(Create(graph), Write(graph_label), run_time=2)
        self.wait(0.5)

        # ---------- 阶段5: 标注零点 ----------
        # lambda = 1 处
        dot1 = Dot(axes.c2p(1, 0), color=RED, radius=0.1)
        label1 = MathTex(r"\lambda = 1", font_size=28, color=RED).next_to(
            dot1, DOWN + LEFT, buff=0.2
        )
        # lambda = 3 处
        dot2 = Dot(axes.c2p(3, 0), color=RED, radius=0.1)
        label2 = MathTex(r"\lambda = 3", font_size=28, color=RED).next_to(
            dot2, DOWN + RIGHT, buff=0.2
        )

        # 从图像到零点的虚线
        dashed1 = DashedLine(
            axes.c2p(1, 1), axes.c2p(1, 0), color=RED, stroke_width=2
        )
        dashed2 = DashedLine(
            axes.c2p(3, 1), axes.c2p(3, 0), color=RED, stroke_width=2
        )

        self.play(
            Create(dashed1), Create(dashed2),
            FadeIn(dot1), FadeIn(dot2),
            Write(label1), Write(label2),
            run_time=1.5,
        )
        self.play(
            FadeIn(eigenvalues_text), Create(ev_box),
            run_time=1,
        )
        self.wait(0.5)

        # 总结文字
        summary = Text(
            "Roots of characteristic polynomial = Eigenvalues",
            font_size=26,
            color=GREEN_A,
        ).to_edge(DOWN)
        self.play(FadeIn(summary), run_time=1)
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 3: 不同类型矩阵的特征值在复平面上的分布
# - 实对称矩阵: 特征值全在实轴
# - 旋转矩阵: 特征值在单位圆上
# - 一般矩阵: 特征值散布在复平面
# ============================================================
class EigenvalueSpectrum(Scene):
    def construct(self):
        title = Text("Eigenvalue Spectrum in the Complex Plane",
                      font_size=30).to_edge(UP)
        self.play(Write(title), run_time=1)

        # 复平面坐标系
        axes = Axes(
            x_range=[-3.5, 3.5, 1],
            y_range=[-2.5, 2.5, 1],
            x_length=9,
            y_length=6,
            axis_config={
                "color": GREY_B,
                "stroke_width": 1.5,
                "include_numbers": True,
                "font_size": 22,
            },
            tips=True,
        ).shift(DOWN * 0.2)
        re_label = axes.get_x_axis_label(
            Text("Re", font_size=24, color=GREY_A), edge=RIGHT, direction=RIGHT
        )
        im_label = axes.get_y_axis_label(
            Text("Im", font_size=24, color=GREY_A), edge=UP, direction=UP
        )

        self.play(Create(axes), Write(re_label), Write(im_label), run_time=1.5)
        self.wait(0.5)

        # ===== Case 1: 实对称矩阵 =====
        self._show_case(
            axes,
            matrix_tex=r"A = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}",
            case_label="Real Symmetric Matrix",
            case_label_color=BLUE_C,
            eigenvalues=[(2, 0), (4, 0)],
            ev_labels=[r"\lambda=2", r"\lambda=4"],
            dot_color=BLUE,
            note="Eigenvalues are always real",
            note_color=BLUE_C,
            highlight_type="real_axis",
            axes_ref=axes,
        )

        # ===== Case 2: 旋转矩阵 =====
        # 旋转 60 度: cos60=0.5, sin60=0.866
        self._show_case(
            axes,
            matrix_tex=r"R = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}",
            case_label="Rotation Matrix (theta = 60 deg)",
            case_label_color=GREEN_C,
            eigenvalues=[(0.5, 0.866), (0.5, -0.866)],
            ev_labels=[r"e^{i\theta}", r"e^{-i\theta}"],
            dot_color=GREEN,
            note="Eigenvalues lie on the unit circle",
            note_color=GREEN_C,
            highlight_type="unit_circle",
            axes_ref=axes,
        )

        # ===== Case 3: 一般矩阵 =====
        self._show_case(
            axes,
            matrix_tex=r"B = \begin{pmatrix} 1 & -2 \\ 3 & 0 \end{pmatrix}",
            case_label="General Matrix",
            case_label_color=ORANGE,
            eigenvalues=[(0.5, 2.398), (0.5, -2.398)],
            ev_labels=[r"\lambda_1", r"\lambda_2"],
            dot_color=ORANGE,
            note="Eigenvalues can be anywhere in the complex plane",
            note_color=ORANGE,
            highlight_type=None,
            axes_ref=axes,
        )

        self.wait(1)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)

    def _show_case(self, axes, matrix_tex, case_label, case_label_color,
                   eigenvalues, ev_labels, dot_color, note, note_color,
                   highlight_type, axes_ref):
        """通用方法：展示一种矩阵类型的特征值分布"""
        # 矩阵标注
        mat = MathTex(matrix_tex, font_size=32).to_corner(UL).shift(DOWN * 0.8)
        mat_box = SurroundingRectangle(mat, color=case_label_color, buff=0.12)
        case_text = Text(case_label, font_size=24, color=case_label_color).next_to(
            mat, DOWN, buff=0.3
        )

        # 高亮区域
        highlight = VGroup()
        if highlight_type == "real_axis":
            # 高亮实轴
            hl_line = Line(
                axes.c2p(-3.5, 0), axes.c2p(3.5, 0),
                color=BLUE, stroke_width=4, stroke_opacity=0.5,
            )
            highlight.add(hl_line)
        elif highlight_type == "unit_circle":
            # 单位圆
            circle = Circle(radius=axes.c2p(1, 0)[0] - axes.c2p(0, 0)[0],
                            color=GREEN, stroke_width=2, stroke_opacity=0.5)
            circle.move_to(axes.c2p(0, 0))
            highlight.add(circle)

        # 特征值点
        dots = VGroup()
        labels = VGroup()
        for (re_val, im_val), lbl in zip(eigenvalues, ev_labels):
            d = Dot(axes.c2p(re_val, im_val), color=dot_color, radius=0.12)
            l = MathTex(lbl, font_size=26, color=dot_color).next_to(d, UR, buff=0.15)
            dots.add(d)
            labels.add(l)

        # 说明文字
        note_text = Text(note, font_size=22, color=note_color).to_edge(DOWN)

        # 播放动画
        self.play(Write(mat), Create(mat_box), Write(case_text), run_time=1.2)
        if highlight.submobjects:
            self.play(Create(highlight), run_time=0.8)
        self.play(
            LaggedStart(*[FadeIn(d, scale=1.5) for d in dots], lag_ratio=0.3),
            run_time=1,
        )
        self.play(
            LaggedStart(*[Write(l) for l in labels], lag_ratio=0.3),
            FadeIn(note_text),
            run_time=1,
        )
        self.wait(2)

        # 清除该 case 的元素
        self.play(
            FadeOut(mat), FadeOut(mat_box), FadeOut(case_text),
            FadeOut(dots), FadeOut(labels),
            FadeOut(highlight), FadeOut(note_text),
            run_time=0.8,
        )


# ============================================================
# Scene 4: 线性动力系统 dx/dt = Ax 的相图
# 三种情况: 稳定节点、不稳定节点、鞍点
# ============================================================
class LinearDynamics(Scene):
    def construct(self):
        title = Text("Phase Portraits of Linear Dynamical Systems",
                      font_size=28).to_edge(UP)
        subtitle = MathTex(
            r"\frac{d\vec{x}}{dt} = A\vec{x}",
            font_size=36,
        ).next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle), run_time=1.5)
        self.wait(1)
        self.play(FadeOut(subtitle), run_time=0.5)

        # ===== Case 1: 稳定节点 (两个负特征值) =====
        self._show_phase_portrait(
            title_text="Stable Node",
            matrix=np.array([[-2, 0.5], [0.5, -1]]),
            matrix_tex=r"A = \begin{pmatrix} -2 & 0.5 \\ 0.5 & -1 \end{pmatrix}",
            eigenvalues_tex=r"\lambda_1 = -2.28, \quad \lambda_2 = -0.72",
            eigenvectors=None,  # 自动计算
            case_color=BLUE_C,
            stream_color_fn=lambda p: BLUE_C,
        )

        # ===== Case 2: 不稳定节点 (两个正特征值) =====
        self._show_phase_portrait(
            title_text="Unstable Node",
            matrix=np.array([[2, 0.5], [0.5, 1]]),
            matrix_tex=r"A = \begin{pmatrix} 2 & 0.5 \\ 0.5 & 1 \end{pmatrix}",
            eigenvalues_tex=r"\lambda_1 = 2.28, \quad \lambda_2 = 0.72",
            eigenvectors=None,
            case_color=RED_C,
            stream_color_fn=lambda p: RED_C,
        )

        # ===== Case 3: 鞍点 (一正一负特征值) =====
        self._show_phase_portrait(
            title_text="Saddle Point",
            matrix=np.array([[1, 2], [2, -1]]),
            matrix_tex=r"A = \begin{pmatrix} 1 & 2 \\ 2 & -1 \end{pmatrix}",
            eigenvalues_tex=r"\lambda_1 = -\sqrt{5} \approx -2.24, \quad \lambda_2 = \sqrt{5} \approx 2.24",
            eigenvectors=None,
            case_color=ORANGE,
            stream_color_fn=lambda p: ORANGE,
        )

        self.wait(1)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)

    def _show_phase_portrait(self, title_text, matrix, matrix_tex,
                              eigenvalues_tex, eigenvectors, case_color,
                              stream_color_fn):
        """绘制一种线性动力系统的相图"""
        A = matrix

        # 计算特征值和特征向量
        evals, evecs = np.linalg.eig(A)

        # 坐标系
        plane = NumberPlane(
            x_range=[-4, 4, 1],
            y_range=[-3, 3, 1],
            x_length=8,
            y_length=5.5,
            background_line_style={
                "stroke_color": GREY_E,
                "stroke_width": 0.8,
                "stroke_opacity": 0.3,
            },
        ).shift(DOWN * 0.3)

        # 标题和矩阵信息
        case_title = Text(title_text, font_size=28, color=case_color).to_edge(UP)
        mat_label = MathTex(matrix_tex, font_size=28).to_corner(UR).shift(
            DOWN * 0.6 + LEFT * 0.2
        )
        ev_label = MathTex(eigenvalues_tex, font_size=22).next_to(
            mat_label, DOWN, buff=0.25
        )
        mat_bg = SurroundingRectangle(
            VGroup(mat_label, ev_label), color=case_color, buff=0.12,
            fill_opacity=0.05,
        )

        # 流线 (使用 ArrowVectorField)
        def vector_field_func(pos):
            x, y = pos[0], pos[1]
            # 将屏幕坐标转换为数值坐标
            coords = plane.p2c(pos)
            v = A @ np.array(coords)
            # 缩放以便可视化，限制长度
            speed = np.sqrt(v[0]**2 + v[1]**2)
            scale = 0.4 / max(speed, 0.5)
            return np.array([v[0] * scale * plane.get_x_unit_size(),
                             v[1] * scale * plane.get_y_unit_size(), 0])

        # 用 StreamLines 代替 ArrowVectorField 来展示流线
        def stream_func(pos):
            coords = plane.p2c(pos)
            v = A @ np.array(coords)
            return np.array([v[0] * plane.get_x_unit_size(),
                             v[1] * plane.get_y_unit_size(), 0])

        stream = StreamLines(
            stream_func,
            x_range=[plane.c2p(-3.5, 0)[0], plane.c2p(3.5, 0)[0], 0.5],
            y_range=[plane.c2p(0, -2.5)[1], plane.c2p(0, 2.5)[1], 0.5],
            stroke_width=1.8,
            max_anchors_per_line=30,
            padding=0.5,
            color=case_color,
            opacity=0.7,
        )

        # 特征向量方向箭头
        eigen_arrows = VGroup()
        eigen_labels_group = VGroup()
        ev_colors = [YELLOW, GREEN_B]
        for i in range(2):
            v = evecs[:, i]
            v_norm = v / np.linalg.norm(v)
            # 正方向
            end_pos = plane.c2p(*(v_norm * 2.8))
            start_pos = plane.c2p(*(- v_norm * 2.8))
            arrow = DashedLine(
                start_pos, end_pos,
                color=ev_colors[i],
                stroke_width=3,
                dash_length=0.15,
            )
            eigen_arrows.add(arrow)

            # 标签
            lbl = MathTex(
                rf"\vec{{v}}_{i+1}",
                font_size=26,
                color=ev_colors[i],
            ).next_to(
                plane.c2p(*(v_norm * 2.5)),
                UR if v_norm[1] >= 0 else DR,
                buff=0.15,
            )
            eigen_labels_group.add(lbl)

        # 播放动画
        self.play(
            Write(case_title),
            Create(plane),
            run_time=1.2,
        )
        self.play(
            Write(mat_label), Write(ev_label), Create(mat_bg),
            run_time=1,
        )
        # 流线动画
        self.play(Create(stream), run_time=2)
        self.wait(0.5)

        # 特征向量方向
        self.play(
            *[Create(a) for a in eigen_arrows],
            *[Write(l) for l in eigen_labels_group],
            run_time=1.2,
        )

        ev_direction_note = Text(
            "Dashed lines = eigenvector directions",
            font_size=20,
            color=GREY_A,
        ).to_edge(DOWN)
        self.play(FadeIn(ev_direction_note), run_time=0.5)
        self.wait(2)

        # 清除
        self.play(
            *[FadeOut(m) for m in [
                case_title, plane, mat_label, ev_label, mat_bg,
                stream, ev_direction_note,
            ]],
            *[FadeOut(a) for a in eigen_arrows],
            *[FadeOut(l) for l in eigen_labels_group],
            run_time=1,
        )
