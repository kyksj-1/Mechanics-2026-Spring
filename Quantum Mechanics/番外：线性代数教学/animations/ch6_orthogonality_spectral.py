"""
第6章：正交性、投影与谱定理 - Manim 可视化动画
==============================================

文件用途：为线性代数教学第6章制作4个可视化场景
  Scene 1 - OrthogonalProjection2D:  二维正交投影
  Scene 2 - GramSchmidtProcess:      Gram-Schmidt 正交化
  Scene 3 - SpectralDecomposition:   对称矩阵的谱分解
  Scene 4 - ProjectionOntoSubspace3D: 三维子空间投影

输出格式：MP4 动画
运行命令（低质量预览）：
  conda activate teaching_env && manim -ql --no_latex_cleanup ch6_orthogonality_spectral.py OrthogonalProjection2D
  conda activate teaching_env && manim -ql --no_latex_cleanup ch6_orthogonality_spectral.py GramSchmidtProcess
  conda activate teaching_env && manim -ql --no_latex_cleanup ch6_orthogonality_spectral.py SpectralDecomposition
  conda activate teaching_env && manim -ql --no_latex_cleanup ch6_orthogonality_spectral.py ProjectionOntoSubspace3D

运行命令（中等质量）：
  conda activate teaching_env && manim -qm --no_latex_cleanup ch6_orthogonality_spectral.py <SceneName>

作者：kyksj-1
"""

from manim import *
import numpy as np


# ============================================================
# 辅助函数
# ============================================================

def create_right_angle_mark(vertex, dir1, dir2, size=0.25, color=WHITE):
    """在 vertex 处创建直角标记，dir1 和 dir2 是两个方向向量（单位向量）"""
    p1 = vertex + dir1 * size
    p2 = vertex + dir2 * size
    p3 = vertex + dir1 * size + dir2 * size
    mark = VGroup(
        Line(p1, p3, color=color, stroke_width=2),
        Line(p2, p3, color=color, stroke_width=2),
    )
    return mark


# ============================================================
# Scene 1: 二维正交投影
# 向量 b = (3, 2) 投影到方向 a = (2, 1) 上
# 投影公式: p = (a^T b / a^T a) * a
# 误差向量: e = b - p 垂直于 a
# ============================================================
class OrthogonalProjection2D(Scene):
    def construct(self):
        # ---------- 阶段1: 搭建坐标平面 ----------
        plane = NumberPlane(
            x_range=[-1, 5, 1],
            y_range=[-1, 4, 1],
            x_length=8,
            y_length=6,
            background_line_style={
                "stroke_color": BLUE_E,
                "stroke_width": 1,
                "stroke_opacity": 0.3,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 2},
        )
        plane.shift(LEFT * 0.5 + DOWN * 0.3)

        # 标题
        title = Text("Orthogonal Projection in 2D", font_size=32).to_edge(UP)

        self.play(Create(plane), Write(title), run_time=1.5)
        self.wait(0.5)

        # ---------- 阶段2: 定义向量并显示 ----------
        # 原始向量
        a_vec = np.array([2, 1])  # 方向向量
        b_vec = np.array([3, 2])  # 被投影的向量

        # 计算投影
        proj_scalar = np.dot(a_vec, b_vec) / np.dot(a_vec, a_vec)  # = 8/5
        p_vec = proj_scalar * a_vec  # 投影向量 p = (16/5, 8/5)
        e_vec = b_vec - p_vec  # 误差向量 e = b - p

        # 绘制方向 a 的延伸线（表示子空间）
        a_unit = a_vec / np.linalg.norm(a_vec)
        line_start = plane.c2p(*(a_unit * (-0.5)))
        line_end = plane.c2p(*(a_unit * 4.5))
        subspace_line = DashedLine(
            line_start, line_end,
            color=BLUE_C,
            stroke_width=2,
            stroke_opacity=0.6,
            dash_length=0.12,
        )
        subspace_label = Text(
            "Subspace spanned by a", font_size=18, color=BLUE_C
        ).next_to(plane.c2p(4, 2), UP, buff=0.15)

        # 向量 a（蓝色）
        arrow_a = Arrow(
            plane.c2p(0, 0), plane.c2p(*a_vec),
            buff=0, stroke_width=4, color=BLUE_C,
            max_tip_length_to_length_ratio=0.15,
        )
        label_a = MathTex(
            r"\vec{a} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}",
            font_size=28, color=BLUE_C,
        ).next_to(arrow_a.get_end(), DR, buff=0.15)

        # 向量 b（黄色）
        arrow_b = Arrow(
            plane.c2p(0, 0), plane.c2p(*b_vec),
            buff=0, stroke_width=4, color=YELLOW,
            max_tip_length_to_length_ratio=0.15,
        )
        label_b = MathTex(
            r"\vec{b} = \begin{pmatrix} 3 \\ 2 \end{pmatrix}",
            font_size=28, color=YELLOW,
        ).next_to(arrow_b.get_end(), UL, buff=0.15)

        # 显示向量 a 和 b
        self.play(
            Create(subspace_line),
            Write(subspace_label),
            run_time=1,
        )
        self.play(
            GrowArrow(arrow_a), Write(label_a),
            run_time=1.2,
        )
        self.play(
            GrowArrow(arrow_b), Write(label_b),
            run_time=1.2,
        )
        self.wait(1)

        # ---------- 阶段3: 展示投影公式 ----------
        formula_box = VGroup(
            MathTex(
                r"\text{proj}_{\vec{a}} \vec{b} = \frac{\vec{a}^T \vec{b}}{\vec{a}^T \vec{a}} \vec{a}",
                font_size=30,
            ),
            MathTex(
                r"= \frac{8}{5} \begin{pmatrix} 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 3.2 \\ 1.6 \end{pmatrix}",
                font_size=28,
            ),
        ).arrange(DOWN, buff=0.2).to_corner(UR).shift(DOWN * 0.8 + LEFT * 0.2)
        formula_bg = SurroundingRectangle(
            formula_box, color=GREEN_C, buff=0.15, fill_opacity=0.08
        )

        self.play(Write(formula_box), Create(formula_bg), run_time=1.5)
        self.wait(1)

        # ---------- 阶段4: 动画投影过程 ----------
        # 投影向量 p（绿色）
        arrow_p = Arrow(
            plane.c2p(0, 0), plane.c2p(*p_vec),
            buff=0, stroke_width=4, color=GREEN_C,
            max_tip_length_to_length_ratio=0.15,
        )
        label_p = MathTex(
            r"\vec{p}", font_size=28, color=GREEN_C,
        ).next_to(arrow_p.get_end(), DOWN, buff=0.2)

        # 从 b 到 p 的投影动画：虚线从 b 末端垂直落到子空间线上
        drop_line = DashedLine(
            plane.c2p(*b_vec), plane.c2p(*p_vec),
            color=GREY_A, stroke_width=2, dash_length=0.1,
        )

        # 先画虚线引导视线
        self.play(Create(drop_line), run_time=1)
        self.play(GrowArrow(arrow_p), Write(label_p), run_time=1.5)
        self.wait(0.8)

        # ---------- 阶段5: 画出误差向量 e ----------
        # 误差向量 e = b - p（红色，从 p 末端到 b 末端）
        arrow_e = Arrow(
            plane.c2p(*p_vec), plane.c2p(*b_vec),
            buff=0, stroke_width=3.5, color=RED_B,
            max_tip_length_to_length_ratio=0.18,
        )
        label_e = MathTex(
            r"\vec{e} = \vec{b} - \vec{p}",
            font_size=26, color=RED_B,
        ).next_to(arrow_e.get_center(), RIGHT, buff=0.2)

        self.play(
            FadeOut(drop_line),
            GrowArrow(arrow_e), Write(label_e),
            run_time=1.5,
        )
        self.wait(0.8)

        # ---------- 阶段6: 直角标记 ----------
        # 在 p 末端标记 e 垂直于 a
        # 方向1: 沿 a 的方向
        a_unit_screen = (
            np.array(plane.c2p(*a_unit)) - np.array(plane.c2p(0, 0))
        )
        a_unit_screen = a_unit_screen / np.linalg.norm(a_unit_screen)
        # 方向2: 沿 e 的方向
        e_unit = e_vec / np.linalg.norm(e_vec)
        e_unit_screen = (
            np.array(plane.c2p(*e_unit)) - np.array(plane.c2p(0, 0))
        )
        e_unit_screen = e_unit_screen / np.linalg.norm(e_unit_screen)

        right_angle = create_right_angle_mark(
            vertex=np.array(plane.c2p(*p_vec)),
            dir1=-a_unit_screen,  # 沿 a 的负方向（回退一点）
            dir2=e_unit_screen,   # 沿 e 的方向
            size=0.2,
            color=WHITE,
        )

        perp_label = MathTex(
            r"\vec{e} \perp \vec{a}",
            font_size=26, color=WHITE,
        ).next_to(right_angle, LEFT, buff=0.3)

        self.play(Create(right_angle), Write(perp_label), run_time=1)
        self.wait(1)

        # ---------- 阶段7: 关键性质总结 ----------
        key_property = VGroup(
            MathTex(
                r"\vec{a}^T (\vec{b} - \vec{p}) = 0",
                font_size=30, color=GREEN_A,
            ),
            Text(
                "The error is orthogonal to the subspace",
                font_size=22, color=GREY_A,
            ),
        ).arrange(DOWN, buff=0.2).to_edge(DOWN)

        self.play(FadeIn(key_property), run_time=1)
        self.wait(2)

        # 淡出所有
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 2: Gram-Schmidt 正交化过程 (2D)
# 输入: v1 = (2, 1), v2 = (1, 2)
# Step 1: u1 = v1
# Step 2: u2 = v2 - proj_{u1}(v2)
# 最终得到正交基 {u1, u2}
# ============================================================
class GramSchmidtProcess(Scene):
    def construct(self):
        # ---------- 阶段1: 坐标系和标题 ----------
        plane = NumberPlane(
            x_range=[-2, 4, 1],
            y_range=[-2, 4, 1],
            x_length=7,
            y_length=7,
            background_line_style={
                "stroke_color": BLUE_E,
                "stroke_width": 1,
                "stroke_opacity": 0.3,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 2},
        )
        plane.shift(LEFT * 1 + DOWN * 0.3)

        title = Text(
            "Gram-Schmidt Orthogonalization", font_size=30
        ).to_edge(UP)

        self.play(Create(plane), Write(title), run_time=1.5)
        self.wait(0.5)

        # ---------- 阶段2: 显示原始向量 v1, v2 ----------
        v1 = np.array([2.0, 1.0])
        v2 = np.array([1.0, 2.0])

        arrow_v1 = Arrow(
            plane.c2p(0, 0), plane.c2p(*v1),
            buff=0, stroke_width=4, color=BLUE_C,
            max_tip_length_to_length_ratio=0.15,
        )
        label_v1 = MathTex(
            r"\vec{v}_1 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}",
            font_size=26, color=BLUE_C,
        ).next_to(arrow_v1.get_end(), DR, buff=0.15)

        arrow_v2 = Arrow(
            plane.c2p(0, 0), plane.c2p(*v2),
            buff=0, stroke_width=4, color=ORANGE,
            max_tip_length_to_length_ratio=0.15,
        )
        label_v2 = MathTex(
            r"\vec{v}_2 = \begin{pmatrix} 1 \\ 2 \end{pmatrix}",
            font_size=26, color=ORANGE,
        ).next_to(arrow_v2.get_end(), UL, buff=0.15)

        # 显示原始向量
        input_label = Text(
            "Input: two non-orthogonal vectors", font_size=22, color=GREY_A,
        ).to_edge(DOWN)

        self.play(
            GrowArrow(arrow_v1), Write(label_v1),
            GrowArrow(arrow_v2), Write(label_v2),
            FadeIn(input_label),
            run_time=1.5,
        )
        self.wait(1.5)

        # ---------- 阶段3: Step 1 - u1 = v1 ----------
        step1_title = Text(
            "Step 1", font_size=24, color=GREEN_C,
        ).to_corner(UR).shift(DOWN * 0.6 + LEFT * 0.3)
        step1_formula = MathTex(
            r"\vec{u}_1 = \vec{v}_1",
            font_size=28,
        ).next_to(step1_title, DOWN, buff=0.2)
        step1_box = SurroundingRectangle(
            VGroup(step1_title, step1_formula),
            color=GREEN_C, buff=0.15, fill_opacity=0.05,
        )

        # 高亮 v1 变为 u1
        arrow_u1 = Arrow(
            plane.c2p(0, 0), plane.c2p(*v1),
            buff=0, stroke_width=5, color=GREEN_C,
            max_tip_length_to_length_ratio=0.15,
        )
        label_u1 = MathTex(
            r"\vec{u}_1", font_size=28, color=GREEN_C,
        ).next_to(arrow_u1.get_end(), DR, buff=0.15)

        self.play(
            FadeOut(input_label),
            Write(step1_title), Write(step1_formula), Create(step1_box),
            run_time=1,
        )
        self.play(
            Transform(arrow_v1, arrow_u1),
            Transform(label_v1, label_u1),
            run_time=1.5,
        )
        self.wait(1)

        # ---------- 阶段4: Step 2 - 计算投影并减去 ----------
        step2_title = Text(
            "Step 2", font_size=24, color=RED_B,
        ).next_to(step1_box, DOWN, buff=0.4)
        step2_formula = VGroup(
            MathTex(
                r"\vec{u}_2 = \vec{v}_2 - \frac{\vec{u}_1^T \vec{v}_2}{\vec{u}_1^T \vec{u}_1} \vec{u}_1",
                font_size=26,
            ),
            MathTex(
                r"= \vec{v}_2 - \frac{4}{5} \vec{u}_1",
                font_size=26,
            ),
        ).arrange(DOWN, buff=0.15).next_to(step2_title, DOWN, buff=0.2)
        step2_box = SurroundingRectangle(
            VGroup(step2_title, step2_formula),
            color=RED_B, buff=0.15, fill_opacity=0.05,
        )

        self.play(
            Write(step2_title), Write(step2_formula), Create(step2_box),
            run_time=1.5,
        )
        self.wait(1)

        # 计算投影分量
        proj_scalar = np.dot(v1, v2) / np.dot(v1, v1)  # = 4/5
        proj_component = proj_scalar * v1  # 投影分量
        u2 = v2 - proj_component  # 正交化后的向量

        # 显示投影分量向量（红色，从原点出发）
        arrow_proj = Arrow(
            plane.c2p(0, 0), plane.c2p(*proj_component),
            buff=0, stroke_width=3, color=RED_B,
            max_tip_length_to_length_ratio=0.18,
        )
        label_proj = MathTex(
            r"\text{proj}_{\vec{u}_1} \vec{v}_2",
            font_size=22, color=RED_B,
        ).next_to(arrow_proj.get_center(), DOWN, buff=0.2)

        # 显示从 v2 末端到投影末端的连线（分解可视化）
        decomp_line = DashedLine(
            plane.c2p(*v2), plane.c2p(*proj_component),
            color=GREY_A, stroke_width=1.5, dash_length=0.1,
        )

        self.play(
            GrowArrow(arrow_proj), Write(label_proj),
            Create(decomp_line),
            run_time=1.5,
        )
        self.wait(1)

        # 动画展示"减去"投影分量的过程：
        # v2 的箭头变形为 u2，同时投影分量淡出
        arrow_u2 = Arrow(
            plane.c2p(0, 0), plane.c2p(*u2),
            buff=0, stroke_width=5, color=TEAL_C,
            max_tip_length_to_length_ratio=0.18,
        )
        label_u2 = MathTex(
            r"\vec{u}_2", font_size=28, color=TEAL_C,
        ).next_to(arrow_u2.get_end(), UL, buff=0.15)

        # 减去过程：先让 v2 闪烁，然后变形为 u2
        self.play(Indicate(arrow_v2, color=ORANGE), run_time=0.8)

        subtract_text = Text(
            "Subtract projection", font_size=22, color=RED_B,
        ).to_edge(DOWN)
        self.play(FadeIn(subtract_text), run_time=0.5)

        self.play(
            Transform(arrow_v2, arrow_u2),
            Transform(label_v2, label_u2),
            FadeOut(arrow_proj), FadeOut(label_proj),
            FadeOut(decomp_line),
            run_time=2,
        )
        self.wait(1)

        # ---------- 阶段5: 验证正交性 ----------
        self.play(FadeOut(subtract_text), run_time=0.3)

        # 直角标记
        u1_unit = v1 / np.linalg.norm(v1)
        u2_unit = u2 / np.linalg.norm(u2)

        u1_screen = np.array(plane.c2p(*u1_unit)) - np.array(plane.c2p(0, 0))
        u1_screen = u1_screen / np.linalg.norm(u1_screen)
        u2_screen = np.array(plane.c2p(*u2_unit)) - np.array(plane.c2p(0, 0))
        u2_screen = u2_screen / np.linalg.norm(u2_screen)

        right_angle = create_right_angle_mark(
            vertex=np.array(plane.c2p(0, 0)),
            dir1=u1_screen,
            dir2=u2_screen,
            size=0.3,
            color=WHITE,
        )

        ortho_check = VGroup(
            MathTex(
                r"\vec{u}_1^T \vec{u}_2 = 0 \;\; \checkmark",
                font_size=28, color=GREEN_A,
            ),
            Text(
                "Orthogonal basis obtained!",
                font_size=22, color=GREEN_A,
            ),
        ).arrange(DOWN, buff=0.2).to_edge(DOWN)

        self.play(Create(right_angle), run_time=0.8)
        self.play(FadeIn(ortho_check), run_time=1)
        self.wait(2)

        # ---------- 阶段6: 展示归一化（可选） ----------
        # 显示归一化后的正交归一基
        u1_norm = v1 / np.linalg.norm(v1)
        u2_norm = u2 / np.linalg.norm(u2)

        normalize_text = Text(
            "Normalize to get orthonormal basis:",
            font_size=22, color=GREY_A,
        ).to_edge(DOWN).shift(UP * 0.8)

        norm_formulas = MathTex(
            r"\hat{e}_1 = \frac{\vec{u}_1}{\|\vec{u}_1\|}, \quad"
            r"\hat{e}_2 = \frac{\vec{u}_2}{\|\vec{u}_2\|}",
            font_size=26,
        ).next_to(normalize_text, DOWN, buff=0.2)

        # 归一化向量箭头
        arrow_e1_norm = Arrow(
            plane.c2p(0, 0), plane.c2p(*u1_norm),
            buff=0, stroke_width=5, color=GREEN_C,
            max_tip_length_to_length_ratio=0.2,
        )
        arrow_e2_norm = Arrow(
            plane.c2p(0, 0), plane.c2p(*u2_norm),
            buff=0, stroke_width=5, color=TEAL_C,
            max_tip_length_to_length_ratio=0.2,
        )

        self.play(FadeOut(ortho_check), run_time=0.5)
        self.play(
            FadeIn(normalize_text), Write(norm_formulas),
            run_time=1,
        )
        self.play(
            Transform(arrow_v1, arrow_e1_norm),
            Transform(arrow_v2, arrow_e2_norm),
            FadeOut(label_v1), FadeOut(label_v2),
            FadeOut(right_angle),
            run_time=1.5,
        )

        # 新的直角标记（归一化后）
        right_angle_2 = create_right_angle_mark(
            vertex=np.array(plane.c2p(0, 0)),
            dir1=u1_screen,
            dir2=u2_screen,
            size=0.25,
            color=YELLOW,
        )
        self.play(Create(right_angle_2), run_time=0.5)
        self.wait(2)

        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 3: 对称矩阵的谱分解
# A = [[2, 1], [1, 2]]
# 特征值: lambda1 = 3, lambda2 = 1
# 特征向量: q1 = (1,1)/sqrt(2), q2 = (1,-1)/sqrt(2)
# 谱分解: A = 3 * q1 q1^T + 1 * q2 q2^T
# ============================================================
class SpectralDecomposition(Scene):
    def construct(self):
        # ---------- 阶段1: 标题和矩阵 ----------
        title = Text(
            "Spectral Decomposition of Symmetric Matrix", font_size=28
        ).to_edge(UP)

        matrix_tex = MathTex(
            r"A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}",
            font_size=36,
        ).to_corner(UL).shift(DOWN * 0.7 + RIGHT * 0.2)
        matrix_bg = SurroundingRectangle(
            matrix_tex, color=BLUE_C, buff=0.12, fill_opacity=0.08
        )

        self.play(Write(title), run_time=1)
        self.play(Write(matrix_tex), Create(matrix_bg), run_time=1)
        self.wait(0.5)

        # ---------- 阶段2: 显示特征分解信息 ----------
        eigen_info = VGroup(
            MathTex(r"\lambda_1 = 3, \quad \vec{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}", font_size=26),
            MathTex(r"\lambda_2 = 1, \quad \vec{q}_2 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix}", font_size=26),
        ).arrange(DOWN, buff=0.2).next_to(matrix_tex, DOWN, buff=0.5)

        self.play(Write(eigen_info), run_time=1.5)
        self.wait(1)

        # ---------- 阶段3: 单位圆和椭圆变换 ----------
        # 清理左上角信息到更紧凑的位置
        self.play(
            FadeOut(title),
            matrix_tex.animate.scale(0.8).to_corner(UL).shift(DOWN * 0.3),
            matrix_bg.animate.scale(0.8).to_corner(UL).shift(DOWN * 0.3),
            eigen_info.animate.scale(0.8).to_corner(UL).shift(DOWN * 1.5),
            run_time=1,
        )

        # 坐标平面
        plane = NumberPlane(
            x_range=[-4.5, 4.5, 1],
            y_range=[-3.5, 3.5, 1],
            x_length=8,
            y_length=6,
            background_line_style={
                "stroke_color": BLUE_E,
                "stroke_width": 0.8,
                "stroke_opacity": 0.25,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1.5},
        )
        plane.shift(RIGHT * 0.8)

        self.play(Create(plane), run_time=1)

        # 矩阵 A
        A = np.array([[2.0, 1.0], [1.0, 2.0]])
        # 特征向量
        q1 = np.array([1, 1]) / np.sqrt(2)
        q2 = np.array([1, -1]) / np.sqrt(2)

        # 画单位圆
        num_pts = 100
        theta_vals = np.linspace(0, 2 * np.pi, num_pts, endpoint=True)
        circle_points = [
            plane.c2p(np.cos(t), np.sin(t)) for t in theta_vals
        ]
        unit_circle = VMobject(color=WHITE, stroke_width=2.5)
        unit_circle.set_points_smoothly(circle_points)

        circle_label = Text(
            "Unit circle", font_size=20, color=WHITE,
        ).next_to(plane.c2p(0, 1.2), UR, buff=0.1)

        self.play(Create(unit_circle), FadeIn(circle_label), run_time=1.5)
        self.wait(1)

        # 变换后的椭圆
        ellipse_points = []
        for t in theta_vals:
            v = np.array([np.cos(t), np.sin(t)])
            Av = A @ v
            ellipse_points.append(plane.c2p(Av[0], Av[1]))
        ellipse = VMobject(color=YELLOW, stroke_width=3)
        ellipse.set_points_smoothly(ellipse_points)

        ellipse_label = Text(
            "Image under A: Av", font_size=20, color=YELLOW,
        ).next_to(plane.c2p(0, 3.2), UP, buff=0.1)

        # 动画：圆变椭圆
        transform_label = Text(
            "Apply A: unit circle -> ellipse",
            font_size=22, color=GREY_A,
        ).to_edge(DOWN)

        self.play(FadeIn(transform_label), run_time=0.5)
        self.play(
            Transform(unit_circle, ellipse),
            Transform(circle_label, ellipse_label),
            run_time=2.5,
        )
        self.wait(1)

        # ---------- 阶段4: 标出特征向量方向（椭圆的主轴） ----------
        self.play(FadeOut(transform_label), run_time=0.3)

        # 特征向量 q1（对应 lambda=3）
        ev1_end = plane.c2p(*(q1 * 3))  # A*q1 = 3*q1, 落在椭圆上
        ev1_neg = plane.c2p(*(q1 * (-3)))
        axis1 = Line(
            ev1_neg, ev1_end,
            color=RED_B, stroke_width=3,
        )
        dot1_pos = Dot(ev1_end, color=RED_B, radius=0.08)
        dot1_neg = Dot(ev1_neg, color=RED_B, radius=0.08)
        axis1_label = MathTex(
            r"\vec{q}_1 \;(\lambda_1 = 3)",
            font_size=22, color=RED_B,
        ).next_to(dot1_pos, UR, buff=0.1)

        # 特征向量 q2（对应 lambda=1）
        ev2_end = plane.c2p(*(q2 * 1))  # A*q2 = 1*q2
        ev2_neg = plane.c2p(*(q2 * (-1)))
        axis2 = Line(
            ev2_neg, ev2_end,
            color=TEAL_C, stroke_width=3,
        )
        dot2_pos = Dot(ev2_end, color=TEAL_C, radius=0.08)
        dot2_neg = Dot(ev2_neg, color=TEAL_C, radius=0.08)
        axis2_label = MathTex(
            r"\vec{q}_2 \;(\lambda_2 = 1)",
            font_size=22, color=TEAL_C,
        ).next_to(dot2_pos, DR, buff=0.1)

        axes_note = Text(
            "Eigenvectors = principal axes of the ellipse",
            font_size=20, color=GREY_A,
        ).to_edge(DOWN)

        self.play(
            Create(axis1), FadeIn(dot1_pos), FadeIn(dot1_neg),
            Write(axis1_label),
            Create(axis2), FadeIn(dot2_pos), FadeIn(dot2_neg),
            Write(axis2_label),
            FadeIn(axes_note),
            run_time=2,
        )
        self.wait(2)

        # ---------- 阶段5: 谱分解公式 ----------
        self.play(
            FadeOut(axes_note),
            FadeOut(unit_circle), FadeOut(circle_label),
            FadeOut(axis1), FadeOut(axis2),
            FadeOut(dot1_pos), FadeOut(dot1_neg),
            FadeOut(dot2_pos), FadeOut(dot2_neg),
            FadeOut(axis1_label), FadeOut(axis2_label),
            FadeOut(plane),
            FadeOut(matrix_tex), FadeOut(matrix_bg),
            FadeOut(eigen_info),
            run_time=1,
        )

        # 谱分解公式居中显示
        spectral_title = Text(
            "Spectral Theorem for Symmetric Matrices",
            font_size=28,
        ).to_edge(UP)

        spectral_formula = MathTex(
            r"A = \sum_{i=1}^{n} \lambda_i \, \vec{q}_i \, \vec{q}_i^T",
            font_size=40,
        ).shift(UP * 1.5)
        spectral_box = SurroundingRectangle(
            spectral_formula, color=YELLOW, buff=0.2
        )

        concrete = MathTex(
            r"A = 3 \cdot \vec{q}_1 \vec{q}_1^T + 1 \cdot \vec{q}_2 \vec{q}_2^T",
            font_size=34,
        ).next_to(spectral_formula, DOWN, buff=0.6)

        self.play(Write(spectral_title), run_time=1)
        self.play(Write(spectral_formula), Create(spectral_box), run_time=1.5)
        self.play(Write(concrete), run_time=1.2)
        self.wait(1)

        # ---------- 阶段6: 展示每个 rank-1 分量 ----------
        # 计算两个 rank-1 矩阵
        # M1 = 3 * q1 @ q1.T = 3 * [[0.5, 0.5], [0.5, 0.5]] = [[1.5, 1.5], [1.5, 1.5]]
        # M2 = 1 * q2 @ q2.T = 1 * [[0.5, -0.5], [-0.5, 0.5]] = [[0.5, -0.5], [-0.5, 0.5]]

        # 在下方并排展示两个分量
        comp1_label = MathTex(
            r"3 \, \vec{q}_1 \vec{q}_1^T = \begin{pmatrix} 1.5 & 1.5 \\ 1.5 & 1.5 \end{pmatrix}",
            font_size=28, color=RED_B,
        )
        comp2_label = MathTex(
            r"1 \, \vec{q}_2 \vec{q}_2^T = \begin{pmatrix} 0.5 & -0.5 \\ -0.5 & 0.5 \end{pmatrix}",
            font_size=28, color=TEAL_C,
        )
        comp_group = VGroup(comp1_label, comp2_label).arrange(
            RIGHT, buff=0.8
        ).next_to(concrete, DOWN, buff=0.6)

        plus_sign = MathTex(r"+", font_size=36).move_to(
            (comp1_label.get_right() + comp2_label.get_left()) / 2
        )

        self.play(
            Write(comp1_label), Write(comp2_label), Write(plus_sign),
            run_time=1.5,
        )
        self.wait(0.5)

        # 下方的小型坐标系展示每个分量的效果
        # 为每个 rank-1 矩阵，显示它将单位圆映射成的图形
        small_plane_left = NumberPlane(
            x_range=[-3, 3, 1], y_range=[-3, 3, 1],
            x_length=3, y_length=3,
            background_line_style={
                "stroke_color": BLUE_E, "stroke_width": 0.5, "stroke_opacity": 0.2
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1},
        ).shift(LEFT * 3 + DOWN * 2.5)

        small_plane_right = NumberPlane(
            x_range=[-3, 3, 1], y_range=[-3, 3, 1],
            x_length=3, y_length=3,
            background_line_style={
                "stroke_color": BLUE_E, "stroke_width": 0.5, "stroke_opacity": 0.2
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1},
        ).shift(RIGHT * 3 + DOWN * 2.5)

        # Rank-1 矩阵将单位圆映射为线段（因为是秩1）
        M1 = 3 * np.outer(q1, q1)
        M2 = 1 * np.outer(q2, q2)

        # M1 将所有点映射到 q1 方向上
        line1_points = []
        for t in theta_vals:
            v = np.array([np.cos(t), np.sin(t)])
            Mv = M1 @ v
            line1_points.append(small_plane_left.c2p(Mv[0], Mv[1]))
        rank1_img_1 = VMobject(color=RED_B, stroke_width=3)
        rank1_img_1.set_points_smoothly(line1_points)

        # M2 将所有点映射到 q2 方向上
        line2_points = []
        for t in theta_vals:
            v = np.array([np.cos(t), np.sin(t)])
            Mv = M2 @ v
            line2_points.append(small_plane_right.c2p(Mv[0], Mv[1]))
        rank1_img_2 = VMobject(color=TEAL_C, stroke_width=3)
        rank1_img_2.set_points_smoothly(line2_points)

        # 小坐标系标签
        label_left = MathTex(
            r"3\,\vec{q}_1\vec{q}_1^T",
            font_size=22, color=RED_B,
        ).next_to(small_plane_left, UP, buff=0.1)
        label_right = MathTex(
            r"1\,\vec{q}_2\vec{q}_2^T",
            font_size=22, color=TEAL_C,
        ).next_to(small_plane_right, UP, buff=0.1)

        self.play(
            Create(small_plane_left), Create(small_plane_right),
            Write(label_left), Write(label_right),
            run_time=1.2,
        )
        self.play(
            Create(rank1_img_1), Create(rank1_img_2),
            run_time=1.5,
        )

        rank1_note = Text(
            "Each rank-1 term projects onto one eigenvector direction",
            font_size=20, color=GREY_A,
        ).to_edge(DOWN)
        self.play(FadeIn(rank1_note), run_time=0.8)
        self.wait(2)

        # ---------- 阶段7: 总结 ----------
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1)

        summary = VGroup(
            Text("Key Insight:", font_size=28, color=YELLOW),
            Text(
                "Every real symmetric matrix can be decomposed",
                font_size=24,
            ),
            Text(
                "into a sum of rank-1 projections along orthogonal eigenvectors.",
                font_size=24,
            ),
            MathTex(
                r"A = Q \Lambda Q^T = \sum_{i} \lambda_i \, \vec{q}_i \, \vec{q}_i^T",
                font_size=36, color=YELLOW,
            ),
        ).arrange(DOWN, buff=0.3).center()
        summary_box = SurroundingRectangle(summary, color=YELLOW, buff=0.3)

        self.play(FadeIn(summary), Create(summary_box), run_time=1.5)
        self.wait(3)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 4: 三维空间中向量到平面的投影
# b = (1, 2, 3)，平面由 a1 = (1, 0, 1)/sqrt(2) 和 a2 = (0, 1, 0) 张成
# 投影公式: p = A(A^T A)^{-1} A^T b
# ============================================================
class ProjectionOntoSubspace3D(ThreeDScene):
    def construct(self):
        # ---------- 阶段1: 标题（平面视角） ----------
        self.set_camera_orientation(phi=0, theta=-PI / 2)

        title = Text(
            "Projection onto a Subspace in 3D", font_size=30
        ).to_edge(UP)
        self.add_fixed_in_frame_mobjects(title)
        self.play(Write(title), run_time=1)
        self.wait(0.5)
        self.play(FadeOut(title), run_time=0.5)
        self.remove(title)

        # ---------- 阶段2: 3D 坐标轴 ----------
        self.set_camera_orientation(
            phi=65 * DEGREES, theta=-50 * DEGREES, zoom=0.8
        )

        axes = ThreeDAxes(
            x_range=[-1, 4, 1],
            y_range=[-1, 4, 1],
            z_range=[-1, 4, 1],
            x_length=5,
            y_length=5,
            z_length=5,
            axis_config={"color": GREY_B, "stroke_width": 1.5},
        )
        x_label = axes.get_x_axis_label(Text("x", font_size=20))
        y_label = axes.get_y_axis_label(Text("y", font_size=20))
        z_label = axes.get_z_axis_label(Text("z", font_size=20))

        self.play(Create(axes), Write(x_label), Write(y_label), Write(z_label), run_time=1.5)
        self.wait(0.5)

        # ---------- 阶段3: 定义平面和向量 ----------
        # 平面的基向量（列向量）
        a1 = np.array([1.0, 0.0, 1.0]) / np.sqrt(2)
        a2 = np.array([0.0, 1.0, 0.0])
        # 被投影的向量
        b = np.array([1.0, 2.0, 3.0])

        # 计算投影
        # A_mat = [a1 | a2]，列向量组成的 3x2 矩阵
        A_mat = np.column_stack([a1, a2])
        # p = A (A^T A)^{-1} A^T b
        ATA = A_mat.T @ A_mat
        ATA_inv = np.linalg.inv(ATA)
        p = A_mat @ ATA_inv @ A_mat.T @ b
        e = b - p  # 误差向量

        # 绘制平面（用 Surface 或 Polygon）
        # 平面参数化: point = s * a1 + t * a2, s,t in [-2, 4]
        plane_surface = Surface(
            lambda u, v: axes.c2p(
                *(u * a1 + v * a2)
            ),
            u_range=[-1.5, 3.5],
            v_range=[-1.5, 3.5],
            resolution=(8, 8),
            fill_opacity=0.2,
            stroke_width=0.5,
            stroke_opacity=0.3,
            checkerboard_colors=[BLUE_C, BLUE_D],
        )

        # 平面基向量箭头
        arrow_a1 = Arrow3D(
            start=axes.c2p(0, 0, 0),
            end=axes.c2p(*(a1 * 2)),
            color=BLUE_C,
        )
        arrow_a2 = Arrow3D(
            start=axes.c2p(0, 0, 0),
            end=axes.c2p(*(a2 * 2)),
            color=BLUE_C,
        )

        # 向量 b
        arrow_b = Arrow3D(
            start=axes.c2p(0, 0, 0),
            end=axes.c2p(*b),
            color=YELLOW,
        )

        # 显示平面和向量
        self.play(Create(plane_surface), run_time=1.5)
        self.play(
            Create(arrow_a1), Create(arrow_a2),
            run_time=1,
        )

        # 固定在画面上的标签
        basis_label = MathTex(
            r"\text{Subspace } W = \text{span}\{\vec{a}_1, \vec{a}_2\}",
            font_size=24, color=BLUE_C,
        ).to_corner(UL).shift(DOWN * 0.5)
        self.add_fixed_in_frame_mobjects(basis_label)
        self.play(Write(basis_label), run_time=0.8)

        self.play(Create(arrow_b), run_time=1)

        b_label = MathTex(
            r"\vec{b} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}",
            font_size=24, color=YELLOW,
        ).to_corner(UR).shift(DOWN * 0.5 + LEFT * 0.2)
        self.add_fixed_in_frame_mobjects(b_label)
        self.play(Write(b_label), run_time=0.8)
        self.wait(1)

        # ---------- 阶段4: 投影动画 ----------
        # 从 b 末端到平面的投影虚线
        drop_line = DashedLine(
            axes.c2p(*b), axes.c2p(*p),
            color=GREY_A, stroke_width=2, dash_length=0.12,
        )

        # 投影向量 p
        arrow_p = Arrow3D(
            start=axes.c2p(0, 0, 0),
            end=axes.c2p(*p),
            color=GREEN_C,
        )

        # 投影点
        p_dot = Dot3D(axes.c2p(*p), color=GREEN_C, radius=0.06)

        self.play(Create(drop_line), run_time=1)
        self.play(Create(arrow_p), FadeIn(p_dot), run_time=1.5)

        p_label = MathTex(
            r"\vec{p} = \text{proj}_W \vec{b}",
            font_size=24, color=GREEN_C,
        ).next_to(b_label, DOWN, buff=0.4)
        self.add_fixed_in_frame_mobjects(p_label)
        self.play(Write(p_label), run_time=0.8)
        self.wait(1)

        # ---------- 阶段5: 误差向量 ----------
        # e = b - p（从 p 到 b 的箭头）
        arrow_e = Arrow3D(
            start=axes.c2p(*p),
            end=axes.c2p(*b),
            color=RED_B,
        )

        self.play(
            FadeOut(drop_line),
            Create(arrow_e),
            run_time=1.5,
        )

        e_label = MathTex(
            r"\vec{e} = \vec{b} - \vec{p} \perp W",
            font_size=24, color=RED_B,
        ).next_to(p_label, DOWN, buff=0.3)
        self.add_fixed_in_frame_mobjects(e_label)
        self.play(Write(e_label), run_time=0.8)
        self.wait(1)

        # ---------- 阶段6: 旋转相机以更好地观察垂直关系 ----------
        # 注意: ThreeDScene 中相机移动需使用 move_camera 方法
        self.move_camera(
            phi=75 * DEGREES,
            theta=-30 * DEGREES,
            run_time=2,
        )
        self.wait(1)

        # ---------- 阶段7: 投影矩阵公式 ----------
        proj_formula = VGroup(
            Text("Projection Matrix:", font_size=22, color=WHITE),
            MathTex(
                r"P = A(A^T A)^{-1} A^T",
                font_size=30, color=GREEN_A,
            ),
            MathTex(
                r"\vec{p} = P \vec{b}",
                font_size=28,
            ),
        ).arrange(DOWN, buff=0.2).to_corner(DL).shift(UP * 0.3 + RIGHT * 0.3)
        proj_box = SurroundingRectangle(
            proj_formula, color=GREEN_A, buff=0.15, fill_opacity=0.1
        )

        self.add_fixed_in_frame_mobjects(proj_formula, proj_box)
        self.play(FadeIn(proj_formula), Create(proj_box), run_time=1.2)
        self.wait(1)

        # ---------- 阶段8: 关键性质 ----------
        key_props = VGroup(
            MathTex(r"P^2 = P", font_size=26, color=YELLOW),
            Text("(idempotent)", font_size=18, color=GREY_A),
            MathTex(r"P^T = P", font_size=26, color=YELLOW),
            Text("(symmetric)", font_size=18, color=GREY_A),
        ).arrange(DOWN, buff=0.15).to_edge(DOWN).shift(UP * 0.2)

        self.add_fixed_in_frame_mobjects(key_props)
        self.play(FadeIn(key_props), run_time=1)

        # 缓慢旋转欣赏
        self.begin_ambient_camera_rotation(rate=0.12)
        self.wait(4)
        self.stop_ambient_camera_rotation()

        # 淡出所有
        all_mobjects = [
            plane_surface, arrow_a1, arrow_a2, arrow_b, arrow_p,
            p_dot, arrow_e, axes, x_label, y_label, z_label,
        ]
        fixed_mobjects = [
            basis_label, b_label, p_label, e_label,
            proj_formula, proj_box, key_props,
        ]
        self.play(
            *[FadeOut(m) for m in all_mobjects],
            run_time=1,
        )
        for m in fixed_mobjects:
            self.remove(m)
