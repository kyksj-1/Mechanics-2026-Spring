"""
第2章：矩阵的相似对角化 - Manim 可视化动画
===========================================

文件用途：为线性代数教学第2章制作4个可视化场景
输出格式：MP4 动画
运行命令：
  conda activate teaching_env && manim -qm ch2_diagonalization_visualizations.py SimilarMatrices
  conda activate teaching_env && manim -qm ch2_diagonalization_visualizations.py DiagonalizationProcess
  conda activate teaching_env && manim -qm ch2_diagonalization_visualizations.py NonDiagonalizable
  conda activate teaching_env && manim -qm ch2_diagonalization_visualizations.py MatrixPowerViaDiagonalization

作者：kyksj-1
"""

from manim import *
import numpy as np


# ============================================================
# 辅助函数
# ============================================================

def apply_matrix_to_point(matrix, point):
    """将2x2矩阵应用于二维坐标点"""
    return matrix @ np.array(point)


def create_basis_arrows(plane, v1, v2, color1=YELLOW, color2=RED_B,
                        label1_tex=None, label2_tex=None, scale=1.5):
    """在平面上创建一对基底向量箭头和标签"""
    arrow1 = Arrow(
        plane.c2p(0, 0), plane.c2p(*(v1 * scale)),
        buff=0, stroke_width=4, color=color1,
        max_tip_length_to_length_ratio=0.15,
    )
    arrow2 = Arrow(
        plane.c2p(0, 0), plane.c2p(*(v2 * scale)),
        buff=0, stroke_width=4, color=color2,
        max_tip_length_to_length_ratio=0.15,
    )
    group = VGroup(arrow1, arrow2)

    labels = VGroup()
    if label1_tex:
        l1 = MathTex(label1_tex, font_size=26, color=color1)
        l1.next_to(arrow1.get_end(), UR, buff=0.1)
        labels.add(l1)
    if label2_tex:
        l2 = MathTex(label2_tex, font_size=26, color=color2)
        l2.next_to(arrow2.get_end(), DR if v2[1] < 0 else UL, buff=0.1)
        labels.add(l2)

    return group, labels


def create_grid_lines(plane, matrix, color=BLUE_C, num_lines=9,
                      x_range=(-4, 4), y_range=(-3, 3), opacity=0.5):
    """创建经过矩阵变换后的网格线"""
    lines = VGroup()
    # 竖直线经过变换
    for i in np.linspace(x_range[0], x_range[1], num_lines):
        pts = []
        for t in np.linspace(y_range[0], y_range[1], 50):
            original = np.array([i, t])
            transformed = matrix @ original
            pts.append(plane.c2p(*transformed))
        line = VMobject(color=color, stroke_width=1.5, stroke_opacity=opacity)
        line.set_points_smoothly(pts)
        lines.add(line)
    # 水平线经过变换
    for j in np.linspace(y_range[0], y_range[1], num_lines):
        pts = []
        for t in np.linspace(x_range[0], x_range[1], 50):
            original = np.array([t, j])
            transformed = matrix @ original
            pts.append(plane.c2p(*transformed))
        line = VMobject(color=color, stroke_width=1.5, stroke_opacity=opacity)
        line.set_points_smoothly(pts)
        lines.add(line)
    return lines


# ============================================================
# Scene 1: 相似矩阵的几何含义
# 同一个线性变换在不同基底下的矩阵表示
# A = [[4,2],[1,3]], 特征值 5,2, 特征向量 (2,1),(-1,1)
# P = [[2,-1],[1,1]], D = [[5,0],[0,2]]
# ============================================================
class SimilarMatrices(Scene):
    def construct(self):
        # ---------- 定义矩阵 ----------
        A = np.array([[4, 2], [1, 3]], dtype=float)
        P = np.array([[2, -1], [1, 1]], dtype=float)
        P_inv = np.linalg.inv(P)
        D = np.diag([5, 2])
        # 特征向量（归一化方向）
        ev1 = np.array([2, 1], dtype=float)
        ev2 = np.array([-1, 1], dtype=float)

        # ===== 阶段1: 引言 =====
        intro_title = Text(
            "Similar Matrices: Same Transformation, Different Bases",
            font_size=28,
        ).to_edge(UP)
        self.play(Write(intro_title), run_time=1.2)

        intro_eq = MathTex(
            r"P^{-1}AP = D",
            font_size=48,
            color=YELLOW,
        ).shift(UP * 0.5)
        intro_desc = Text(
            "A and D represent the same linear transformation\n"
            "in different coordinate systems",
            font_size=22, color=GREY_A,
            line_spacing=1.3,
        ).next_to(intro_eq, DOWN, buff=0.6)

        self.play(Write(intro_eq), run_time=1.5)
        self.play(FadeIn(intro_desc), run_time=1)
        self.wait(2)
        self.play(FadeOut(intro_eq), FadeOut(intro_desc), run_time=0.8)

        # ===== 阶段2: 标准基下的变换 (矩阵 A) =====
        self.play(
            intro_title.animate.become(
                Text("Standard Basis: Matrix A", font_size=28).to_edge(UP)
            ),
            run_time=0.8,
        )

        # 左侧：标准基网格（变换前）
        plane_left = NumberPlane(
            x_range=[-4, 4, 1], y_range=[-3, 3, 1],
            x_length=5.5, y_length=4,
            background_line_style={
                "stroke_color": BLUE_E, "stroke_width": 1, "stroke_opacity": 0.4,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1.5},
        ).shift(LEFT * 3.2 + DOWN * 0.4)

        mat_a_label = MathTex(
            r"A = \begin{pmatrix} 4 & 2 \\ 1 & 3 \end{pmatrix}",
            font_size=30,
        ).to_corner(UL).shift(DOWN * 0.7 + RIGHT * 0.2)
        mat_a_box = SurroundingRectangle(mat_a_label, color=BLUE, buff=0.1,
                                         fill_opacity=0.08)

        # 标准基向量
        e1 = np.array([1, 0], dtype=float)
        e2 = np.array([0, 1], dtype=float)
        std_arrows, std_labels = create_basis_arrows(
            plane_left, e1, e2, TEAL_B, LIGHT_PINK,
            r"\vec{e}_1", r"\vec{e}_2", scale=1.2,
        )

        # 标准基网格
        identity = np.eye(2)
        std_grid = create_grid_lines(
            plane_left, identity, BLUE_D, num_lines=7,
            x_range=(-3, 3), y_range=(-2.5, 2.5), opacity=0.35,
        )

        self.play(Create(plane_left), run_time=1)
        self.play(
            FadeIn(mat_a_label), Create(mat_a_box),
            Create(std_grid),
            run_time=1,
        )
        self.play(
            *[GrowArrow(a) for a in std_arrows],
            *[Write(l) for l in std_labels],
            run_time=1,
        )
        self.wait(0.5)

        # 变换后的网格
        transformed_grid = create_grid_lines(
            plane_left, A, YELLOW_C, num_lines=7,
            x_range=(-3, 3), y_range=(-2.5, 2.5), opacity=0.5,
        )

        # 变换后的基向量
        Ae1 = A @ e1
        Ae2 = A @ e2
        trans_arrows = VGroup(
            Arrow(plane_left.c2p(0, 0), plane_left.c2p(*(Ae1 * 1.2)),
                  buff=0, stroke_width=4, color=TEAL_B,
                  max_tip_length_to_length_ratio=0.12),
            Arrow(plane_left.c2p(0, 0), plane_left.c2p(*(Ae2 * 1.2)),
                  buff=0, stroke_width=4, color=LIGHT_PINK,
                  max_tip_length_to_length_ratio=0.12),
        )

        apply_label = Text("Apply A", font_size=24, color=YELLOW_C).next_to(
            plane_left, DOWN, buff=0.2
        )

        self.play(FadeIn(apply_label), run_time=0.5)
        self.play(
            Transform(std_grid, transformed_grid),
            Transform(std_arrows[0], trans_arrows[0]),
            Transform(std_arrows[1], trans_arrows[1]),
            run_time=2.5,
        )
        self.wait(1)

        note_a = Text(
            "Grid is sheared and stretched",
            font_size=20, color=GREY_A,
        ).next_to(apply_label, DOWN, buff=0.15)
        self.play(FadeIn(note_a), run_time=0.8)
        self.wait(1)

        # ===== 阶段3: 特征基下的变换 (对角矩阵 D) =====
        # 右侧面板
        self.play(
            intro_title.animate.become(
                Text("Eigenbasis: Diagonal Matrix D", font_size=28).to_edge(UP)
            ),
            run_time=0.8,
        )

        plane_right = NumberPlane(
            x_range=[-4, 4, 1], y_range=[-3, 3, 1],
            x_length=5.5, y_length=4,
            background_line_style={
                "stroke_color": GREEN_E, "stroke_width": 1, "stroke_opacity": 0.4,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1.5},
        ).shift(RIGHT * 3.2 + DOWN * 0.4)

        mat_d_label = MathTex(
            r"D = \begin{pmatrix} 5 & 0 \\ 0 & 2 \end{pmatrix}",
            font_size=30,
        ).to_corner(UR).shift(DOWN * 0.7 + LEFT * 0.2)
        mat_d_box = SurroundingRectangle(mat_d_label, color=GREEN, buff=0.1,
                                         fill_opacity=0.08)

        # 在特征基中，坐标轴方向就是特征向量
        # 显示的是特征基坐标下的网格
        ev1_dir = ev1 / np.linalg.norm(ev1)
        ev2_dir = ev2 / np.linalg.norm(ev2)
        eigen_arrows, eigen_labels = create_basis_arrows(
            plane_right, np.array([1, 0]), np.array([0, 1]),
            YELLOW, RED_B,
            r"\vec{v}_1\text{ dir}", r"\vec{v}_2\text{ dir}", scale=1.2,
        )

        eigen_grid = create_grid_lines(
            plane_right, identity, GREEN_D, num_lines=7,
            x_range=(-3, 3), y_range=(-2.5, 2.5), opacity=0.35,
        )

        self.play(Create(plane_right), run_time=1)
        self.play(
            FadeIn(mat_d_label), Create(mat_d_box),
            Create(eigen_grid),
            run_time=1,
        )
        self.play(
            *[GrowArrow(a) for a in eigen_arrows],
            *[Write(l) for l in eigen_labels],
            run_time=1,
        )
        self.wait(0.5)

        # D 的变换：纯拉伸
        transformed_eigen_grid = create_grid_lines(
            plane_right, D, YELLOW_C, num_lines=7,
            x_range=(-3, 3), y_range=(-2.5, 2.5), opacity=0.5,
        )

        De1 = D @ np.array([1, 0])
        De2 = D @ np.array([0, 1])
        diag_trans_arrows = VGroup(
            Arrow(plane_right.c2p(0, 0), plane_right.c2p(*(De1 * 1.2)),
                  buff=0, stroke_width=4, color=YELLOW,
                  max_tip_length_to_length_ratio=0.12),
            Arrow(plane_right.c2p(0, 0), plane_right.c2p(*(De2 * 1.2)),
                  buff=0, stroke_width=4, color=RED_B,
                  max_tip_length_to_length_ratio=0.12),
        )

        apply_label_d = Text("Apply D", font_size=24, color=YELLOW_C).next_to(
            plane_right, DOWN, buff=0.2
        )

        self.play(FadeIn(apply_label_d), run_time=0.5)
        self.play(
            Transform(eigen_grid, transformed_eigen_grid),
            Transform(eigen_arrows[0], diag_trans_arrows[0]),
            Transform(eigen_arrows[1], diag_trans_arrows[1]),
            run_time=2.5,
        )
        self.wait(0.5)

        note_d = Text(
            "Grid is only stretched along axes (pure scaling)",
            font_size=20, color=GREY_A,
        ).next_to(apply_label_d, DOWN, buff=0.15)
        self.play(FadeIn(note_d), run_time=0.8)
        self.wait(1.5)

        # ===== 阶段4: 连接公式 =====
        self.play(
            intro_title.animate.become(
                Text("The Change-of-Basis Relation", font_size=28).to_edge(UP)
            ),
            run_time=0.8,
        )

        # 中央公式
        central_eq = MathTex(
            r"P^{-1}", r"A", r"P", r"=", r"D",
            font_size=44,
        ).shift(DOWN * 2.8)
        central_eq[0].set_color(ORANGE)   # P^{-1}
        central_eq[1].set_color(BLUE)     # A
        central_eq[2].set_color(ORANGE)   # P
        central_eq[4].set_color(GREEN)    # D

        p_label = MathTex(
            r"P = \begin{pmatrix} 2 & -1 \\ 1 & 1 \end{pmatrix}",
            font_size=26, color=ORANGE,
        ).next_to(central_eq, DOWN, buff=0.3)
        p_note = Text(
            "P = [eigenvectors as columns]",
            font_size=18, color=ORANGE,
        ).next_to(p_label, DOWN, buff=0.1)

        self.play(FadeOut(note_a), FadeOut(note_d), run_time=0.5)
        self.play(Write(central_eq), run_time=1.5)
        self.play(Write(p_label), FadeIn(p_note), run_time=1)
        self.wait(1)

        # 连接箭头：从左图到右图
        connect_arrow = CurvedArrow(
            plane_left.get_right() + UP * 0.5,
            plane_right.get_left() + UP * 0.5,
            color=ORANGE, stroke_width=3,
        )
        connect_label = MathTex(
            r"\text{Change basis by } P",
            font_size=22, color=ORANGE,
        ).next_to(connect_arrow, UP, buff=0.1)

        self.play(Create(connect_arrow), Write(connect_label), run_time=1.2)
        self.wait(2)

        # ===== 阶段5: 总结 =====
        summary_box = VGroup(
            Text("Key Insight:", font_size=24, color=WHITE, weight=BOLD),
            Text(
                "Diagonalization = finding the basis where the\n"
                "transformation becomes pure scaling",
                font_size=20, color=GREY_A,
                line_spacing=1.3,
            ),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        summary_rect = SurroundingRectangle(
            summary_box, color=YELLOW, buff=0.2, fill_opacity=0.1,
        )
        summary_group = VGroup(summary_rect, summary_box).to_edge(DOWN).shift(UP * 0.2)

        self.play(FadeOut(p_label), FadeOut(p_note), FadeOut(central_eq), run_time=0.5)
        self.play(FadeIn(summary_group), run_time=1)
        self.wait(2.5)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 2: 对角化求解过程动画
# 矩阵 A = [[4,2],[1,3]]
# 特征值: 5, 2
# 特征向量: (2,1), (-1,1)
# P = [[2,-1],[1,1]], D = [[5,0],[0,2]]
# ============================================================
class DiagonalizationProcess(Scene):
    def construct(self):
        A = np.array([[4, 2], [1, 3]], dtype=float)

        # ===== Step 1: 展示矩阵 A =====
        title = Text("Diagonalization Step by Step", font_size=30).to_edge(UP)
        self.play(Write(title), run_time=1)

        step1_label = Text("Step 1: Start with matrix A", font_size=22,
                           color=BLUE_C).next_to(title, DOWN, buff=0.3)
        mat_a = MathTex(
            r"A = \begin{pmatrix} 4 & 2 \\ 1 & 3 \end{pmatrix}",
            font_size=44,
        )
        self.play(Write(step1_label), run_time=0.8)
        self.play(Write(mat_a), run_time=1.2)
        self.wait(1)

        # ===== Step 2: 求特征值 =====
        step2_label = Text("Step 2: Find eigenvalues", font_size=22,
                           color=GREEN_C)
        self.play(
            FadeOut(step1_label),
            mat_a.animate.scale(0.8).to_corner(UL).shift(DOWN * 0.7),
            run_time=0.8,
        )
        step2_label.next_to(title, DOWN, buff=0.3)
        self.play(Write(step2_label), run_time=0.8)

        # 特征方程推导
        char_eq_1 = MathTex(
            r"\det(A - \lambda I) = 0",
            font_size=36,
        ).shift(UP * 0.8)

        char_eq_2 = MathTex(
            r"\det \begin{pmatrix} 4 - \lambda & 2 \\ 1 & 3 - \lambda \end{pmatrix} = 0",
            font_size=34,
        ).next_to(char_eq_1, DOWN, buff=0.4)

        char_eq_3 = MathTex(
            r"(4-\lambda)(3-\lambda) - 2 = 0",
            font_size=34,
        ).next_to(char_eq_2, DOWN, buff=0.3)

        char_eq_4 = MathTex(
            r"\lambda^2 - 7\lambda + 10 = 0",
            font_size=34,
        ).next_to(char_eq_3, DOWN, buff=0.3)

        char_eq_5 = MathTex(
            r"(\lambda - 5)(\lambda - 2) = 0",
            font_size=36, color=YELLOW,
        ).next_to(char_eq_4, DOWN, buff=0.3)

        self.play(Write(char_eq_1), run_time=1)
        self.wait(0.5)
        self.play(Write(char_eq_2), run_time=1.2)
        self.wait(0.5)
        self.play(Write(char_eq_3), run_time=1)
        self.play(Write(char_eq_4), run_time=1)
        self.play(Write(char_eq_5), run_time=1)
        self.wait(0.5)

        # 特征值结果
        eigenvalues_result = MathTex(
            r"\lambda_1 = 5, \quad \lambda_2 = 2",
            font_size=38, color=YELLOW,
        ).shift(DOWN * 2.2)
        ev_box = SurroundingRectangle(eigenvalues_result, color=YELLOW, buff=0.15)

        self.play(Write(eigenvalues_result), Create(ev_box), run_time=1)
        self.wait(1.5)

        # 清理
        self.play(
            FadeOut(char_eq_1), FadeOut(char_eq_2), FadeOut(char_eq_3),
            FadeOut(char_eq_4), FadeOut(char_eq_5),
            eigenvalues_result.animate.scale(0.7).next_to(mat_a, DOWN, buff=0.3),
            FadeOut(ev_box),
            run_time=0.8,
        )

        # ===== Step 3: 求特征向量 =====
        step3_label = Text("Step 3: Find eigenvectors", font_size=22,
                           color=ORANGE)
        self.play(FadeOut(step2_label), run_time=0.3)
        step3_label.next_to(title, DOWN, buff=0.3)
        self.play(Write(step3_label), run_time=0.8)

        # lambda1 = 5
        ev1_title = MathTex(r"\lambda_1 = 5:", font_size=30, color=YELLOW
                            ).shift(UP * 1.2 + LEFT * 3)
        ev1_eq = MathTex(
            r"(A - 5I)\vec{v} = \begin{pmatrix} -1 & 2 \\ 1 & -2 \end{pmatrix}"
            r"\vec{v} = \vec{0}",
            font_size=28,
        ).next_to(ev1_title, DOWN, buff=0.3)
        ev1_result = MathTex(
            r"\vec{v}_1 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}",
            font_size=32, color=YELLOW,
        ).next_to(ev1_eq, DOWN, buff=0.3)

        # lambda2 = 2
        ev2_title = MathTex(r"\lambda_2 = 2:", font_size=30, color=RED_B
                            ).shift(UP * 1.2 + RIGHT * 3)
        ev2_eq = MathTex(
            r"(A - 2I)\vec{v} = \begin{pmatrix} 2 & 2 \\ 1 & 1 \end{pmatrix}"
            r"\vec{v} = \vec{0}",
            font_size=28,
        ).next_to(ev2_title, DOWN, buff=0.3)
        ev2_result = MathTex(
            r"\vec{v}_2 = \begin{pmatrix} -1 \\ 1 \end{pmatrix}",
            font_size=32, color=RED_B,
        ).next_to(ev2_eq, DOWN, buff=0.3)

        self.play(Write(ev1_title), Write(ev2_title), run_time=0.8)
        self.play(Write(ev1_eq), Write(ev2_eq), run_time=1.5)
        self.play(Write(ev1_result), Write(ev2_result), run_time=1.2)
        self.wait(1.5)

        # 清理
        self.play(
            FadeOut(ev1_title), FadeOut(ev1_eq),
            FadeOut(ev2_title), FadeOut(ev2_eq),
            ev1_result.animate.scale(0.8).next_to(eigenvalues_result, DOWN, buff=0.2),
            ev2_result.animate.scale(0.8).next_to(eigenvalues_result, DOWN, buff=0.6),
            run_time=0.8,
        )

        # ===== Step 4: 构建 P 和 D =====
        step4_label = Text("Step 4: Construct P and D", font_size=22,
                           color=TEAL_C)
        self.play(FadeOut(step3_label), run_time=0.3)
        step4_label.next_to(title, DOWN, buff=0.3)
        self.play(Write(step4_label), run_time=0.8)

        p_mat = MathTex(
            r"P = [\vec{v}_1 \mid \vec{v}_2] = "
            r"\begin{pmatrix} 2 & -1 \\ 1 & 1 \end{pmatrix}",
            font_size=34,
        ).shift(UP * 0.5)

        d_mat = MathTex(
            r"D = \begin{pmatrix} \lambda_1 & 0 \\ 0 & \lambda_2 \end{pmatrix} = "
            r"\begin{pmatrix} 5 & 0 \\ 0 & 2 \end{pmatrix}",
            font_size=34,
        ).next_to(p_mat, DOWN, buff=0.5)

        self.play(Write(p_mat), run_time=1.5)
        self.play(Write(d_mat), run_time=1.5)
        self.wait(1)

        # 验证公式
        verify_eq = MathTex(
            r"P^{-1}AP = D \quad \Longleftrightarrow \quad AP = PD",
            font_size=34, color=GREEN_C,
        ).next_to(d_mat, DOWN, buff=0.6)
        verify_box = SurroundingRectangle(verify_eq, color=GREEN, buff=0.15)

        self.play(Write(verify_eq), Create(verify_box), run_time=1.5)
        self.wait(1.5)

        # 清理进入几何演示
        self.play(
            *[FadeOut(m) for m in self.mobjects],
            run_time=1,
        )

        # ===== Step 5: 换基的几何含义 =====
        geo_title = Text(
            "Geometric Meaning: Change of Basis",
            font_size=28,
        ).to_edge(UP)
        self.play(Write(geo_title), run_time=1)

        # 数值平面
        plane = NumberPlane(
            x_range=[-5, 5, 1], y_range=[-4, 4, 1],
            background_line_style={
                "stroke_color": BLUE_E, "stroke_width": 1, "stroke_opacity": 0.3,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1.5},
        ).scale(0.75)

        self.play(Create(plane), run_time=1)

        # 显示标准基网格
        identity = np.eye(2)
        std_grid = create_grid_lines(
            plane, identity, BLUE_D, num_lines=9,
            x_range=(-4, 4), y_range=(-3.5, 3.5), opacity=0.3,
        )
        self.play(Create(std_grid), run_time=1)

        # 标准基向量
        e1_arrow = Arrow(
            plane.c2p(0, 0), plane.c2p(1.5, 0),
            buff=0, color=TEAL_B, stroke_width=3,
            max_tip_length_to_length_ratio=0.15,
        )
        e2_arrow = Arrow(
            plane.c2p(0, 0), plane.c2p(0, 1.5),
            buff=0, color=LIGHT_PINK, stroke_width=3,
            max_tip_length_to_length_ratio=0.15,
        )
        e1_label = MathTex(r"\vec{e}_1", font_size=24, color=TEAL_B).next_to(
            e1_arrow.get_end(), DOWN, buff=0.1
        )
        e2_label = MathTex(r"\vec{e}_2", font_size=24, color=LIGHT_PINK).next_to(
            e2_arrow.get_end(), LEFT, buff=0.1
        )

        # 特征基向量
        ev1 = np.array([2, 1], dtype=float)
        ev2 = np.array([-1, 1], dtype=float)
        ev1_norm = ev1 / np.linalg.norm(ev1)
        ev2_norm = ev2 / np.linalg.norm(ev2)

        v1_arrow = Arrow(
            plane.c2p(0, 0), plane.c2p(*(ev1_norm * 1.8)),
            buff=0, color=YELLOW, stroke_width=3,
            max_tip_length_to_length_ratio=0.15,
        )
        v2_arrow = Arrow(
            plane.c2p(0, 0), plane.c2p(*(ev2_norm * 1.8)),
            buff=0, color=RED_B, stroke_width=3,
            max_tip_length_to_length_ratio=0.15,
        )
        v1_label = MathTex(r"\vec{v}_1", font_size=24, color=YELLOW).next_to(
            v1_arrow.get_end(), UR, buff=0.1
        )
        v2_label = MathTex(r"\vec{v}_2", font_size=24, color=RED_B).next_to(
            v2_arrow.get_end(), UL, buff=0.1
        )

        # 先展示标准基
        std_basis_label = Text("Standard basis", font_size=22,
                               color=TEAL_B).to_edge(DOWN)
        self.play(
            GrowArrow(e1_arrow), GrowArrow(e2_arrow),
            Write(e1_label), Write(e2_label),
            FadeIn(std_basis_label),
            run_time=1.2,
        )
        self.wait(1)

        # 过渡到特征基
        eigen_basis_label = Text("Eigenbasis", font_size=22,
                                 color=YELLOW).to_edge(DOWN)
        self.play(
            ReplacementTransform(e1_arrow, v1_arrow),
            ReplacementTransform(e2_arrow, v2_arrow),
            ReplacementTransform(e1_label, v1_label),
            ReplacementTransform(e2_label, v2_label),
            ReplacementTransform(std_basis_label, eigen_basis_label),
            run_time=2,
        )
        self.wait(0.5)

        # 画出特征基构成的斜网格线
        P = np.array([[2, -1], [1, 1]], dtype=float) / np.sqrt(5)  # 归一化列
        P_actual = np.array([[2, -1], [1, 1]], dtype=float)
        # 实际上用原始P的列方向构建网格
        eigen_grid = VGroup()
        for i in np.linspace(-3, 3, 7):
            # 沿 v1 方向的线
            pts = [plane.c2p(*(ev1_norm * t + ev2_norm * i))
                   for t in np.linspace(-4, 4, 30)]
            line = VMobject(color=YELLOW, stroke_width=1, stroke_opacity=0.4)
            line.set_points_smoothly(pts)
            eigen_grid.add(line)
            # 沿 v2 方向的线
            pts2 = [plane.c2p(*(ev2_norm * t + ev1_norm * i))
                    for t in np.linspace(-4, 4, 30)]
            line2 = VMobject(color=RED_B, stroke_width=1, stroke_opacity=0.4)
            line2.set_points_smoothly(pts2)
            eigen_grid.add(line2)

        self.play(Create(eigen_grid), run_time=1.5)
        self.wait(0.5)

        change_note = Text(
            "In this basis, the transformation is just scaling",
            font_size=20, color=GREEN_A,
        ).next_to(eigen_basis_label, UP, buff=0.15)
        self.play(FadeIn(change_note), run_time=0.8)

        # 在特征方向上展示拉伸
        scale_v1 = Arrow(
            plane.c2p(0, 0), plane.c2p(*(ev1_norm * 1.8 * 2.5)),
            buff=0, color=YELLOW, stroke_width=5,
            max_tip_length_to_length_ratio=0.1,
        )
        scale_v2 = Arrow(
            plane.c2p(0, 0), plane.c2p(*(ev2_norm * 1.8 * 1.3)),
            buff=0, color=RED_B, stroke_width=5,
            max_tip_length_to_length_ratio=0.1,
        )
        scale_label1 = MathTex(r"\times 5", font_size=26, color=YELLOW).next_to(
            scale_v1.get_end(), UR, buff=0.1
        )
        scale_label2 = MathTex(r"\times 2", font_size=26, color=RED_B).next_to(
            scale_v2.get_end(), UL, buff=0.1
        )

        self.play(
            Transform(v1_arrow, scale_v1),
            Transform(v2_arrow, scale_v2),
            FadeIn(scale_label1), FadeIn(scale_label2),
            run_time=2,
        )
        self.wait(2)

        # 最终公式
        final_eq = MathTex(
            r"A = PDP^{-1}",
            font_size=42, color=GREEN_C,
        ).to_corner(DR).shift(UP * 0.5 + LEFT * 0.5)
        final_box = SurroundingRectangle(final_eq, color=GREEN, buff=0.15)
        self.play(Write(final_eq), Create(final_box), run_time=1.2)
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 3: 不可对角化矩阵
# A = [[3,1],[0,3]], 特征值 3 (代数重数2, 几何重数1)
# 只有一个特征方向 (1,0)
# ============================================================
class NonDiagonalizable(Scene):
    def construct(self):
        A_defective = np.array([[3, 1], [0, 3]], dtype=float)
        A_good = np.array([[4, 2], [1, 3]], dtype=float)  # 可对角化的对比

        # ===== 阶段1: 引入问题 =====
        title = Text("Non-Diagonalizable Matrix", font_size=30).to_edge(UP)
        self.play(Write(title), run_time=1)

        mat_label = MathTex(
            r"A = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}",
            font_size=40,
        ).shift(UP * 1)
        self.play(Write(mat_label), run_time=1.2)
        self.wait(0.5)

        # 分析特征方程
        char_eq = MathTex(
            r"\det(A - \lambda I) = (3-\lambda)^2 = 0",
            font_size=34,
        ).next_to(mat_label, DOWN, buff=0.5)
        eigenvalue_result = MathTex(
            r"\lambda = 3 \quad \text{(algebraic multiplicity 2)}",
            font_size=30, color=YELLOW,
        ).next_to(char_eq, DOWN, buff=0.3)

        self.play(Write(char_eq), run_time=1.2)
        self.play(Write(eigenvalue_result), run_time=1)
        self.wait(1)

        # 求特征向量
        ev_eq = MathTex(
            r"(A - 3I)\vec{v} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}\vec{v} = \vec{0}",
            font_size=32,
        ).next_to(eigenvalue_result, DOWN, buff=0.4)
        ev_result = MathTex(
            r"\text{Only } \vec{v}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}"
            r"\quad \text{(geometric multiplicity 1)}",
            font_size=28, color=RED_B,
        ).next_to(ev_eq, DOWN, buff=0.3)

        self.play(Write(ev_eq), run_time=1.2)
        self.play(Write(ev_result), run_time=1)
        self.wait(1)

        # 关键判断
        conclusion = MathTex(
            r"\text{Geometric mult.} < \text{Algebraic mult.}",
            r"\implies \text{NOT diagonalizable!}",
            font_size=28,
        ).next_to(ev_result, DOWN, buff=0.4)
        conclusion[1].set_color(RED)
        conclusion_box = SurroundingRectangle(conclusion, color=RED, buff=0.15)

        self.play(Write(conclusion), Create(conclusion_box), run_time=1.5)
        self.wait(2)

        # 清理
        self.play(*[FadeOut(m) for m in self.mobjects if m != title], run_time=0.8)

        # ===== 阶段2: 几何展示 - 不可对角化矩阵的变换效果 =====
        self.play(
            title.animate.become(
                Text("Geometric Effect: Shear + Scale", font_size=28).to_edge(UP)
            ),
            run_time=0.8,
        )

        plane = NumberPlane(
            x_range=[-5, 5, 1], y_range=[-4, 4, 1],
            background_line_style={
                "stroke_color": BLUE_E, "stroke_width": 1, "stroke_opacity": 0.3,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1.5},
        ).scale(0.7)

        self.play(Create(plane), run_time=1)

        # 矩阵标签
        mat_corner = MathTex(
            r"A = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}",
            font_size=28,
        ).to_corner(UL).shift(DOWN * 0.7)
        mat_bg = SurroundingRectangle(mat_corner, color=RED_B, buff=0.1,
                                      fill_opacity=0.08)
        self.play(FadeIn(mat_corner), Create(mat_bg), run_time=0.8)

        # 画一组向量（变换前）
        num_vectors = 12
        angles = np.linspace(0, 2 * np.pi, num_vectors, endpoint=False)
        radius = 2.0
        colors = color_gradient([BLUE_C, GREEN_C, TEAL_C], num_vectors)

        orig_vecs = VGroup()
        trans_vecs = VGroup()
        for i, angle in enumerate(angles):
            orig_pt = np.array([radius * np.cos(angle), radius * np.sin(angle)])
            trans_pt = A_defective @ orig_pt
            ov = Arrow(
                plane.c2p(0, 0), plane.c2p(*orig_pt),
                buff=0, stroke_width=3, color=colors[i],
                max_tip_length_to_length_ratio=0.15,
            )
            tv = Arrow(
                plane.c2p(0, 0), plane.c2p(*trans_pt),
                buff=0, stroke_width=3, color=colors[i],
                max_tip_length_to_length_ratio=0.15,
            )
            orig_vecs.add(ov)
            trans_vecs.add(tv)

        self.play(
            LaggedStart(*[GrowArrow(v) for v in orig_vecs], lag_ratio=0.06),
            run_time=1.5,
        )
        self.wait(0.5)

        # 施加变换
        apply_text = Text("Apply A", font_size=22, color=YELLOW_C).to_edge(DOWN)
        self.play(FadeIn(apply_text), run_time=0.3)
        self.play(
            *[Transform(orig_vecs[i], trans_vecs[i]) for i in range(num_vectors)],
            run_time=2.5,
        )
        self.wait(0.5)

        # 高亮唯一的特征方向
        eigen_dir = Arrow(
            plane.c2p(-3, 0), plane.c2p(3, 0),
            buff=0, stroke_width=5, color=YELLOW,
            max_tip_length_to_length_ratio=0.08,
        )
        eigen_label = MathTex(
            r"\text{Only eigendirection: } (1,0)",
            font_size=24, color=YELLOW,
        ).next_to(eigen_dir, UP, buff=0.3).shift(RIGHT * 1)

        shear_note = Text(
            "Circle -> Ellipse (shear effect from off-diagonal 1)",
            font_size=18, color=GREY_A,
        ).next_to(apply_text, UP, buff=0.15)

        self.play(
            GrowArrow(eigen_dir),
            Write(eigen_label),
            FadeIn(shear_note),
            run_time=1.5,
        )
        self.wait(2)

        # ===== 阶段3: 对比可对角化矩阵 =====
        self.play(
            *[FadeOut(m) for m in self.mobjects if m != title],
            run_time=0.8,
        )
        self.play(
            title.animate.become(
                Text("Comparison: Diagonalizable vs Non-Diagonalizable",
                     font_size=26).to_edge(UP)
            ),
            run_time=0.8,
        )

        # 左右对比布局
        # 左侧: 可对角化
        plane_l = NumberPlane(
            x_range=[-4, 4, 1], y_range=[-3, 3, 1],
            x_length=5, y_length=3.8,
            background_line_style={
                "stroke_color": BLUE_E, "stroke_width": 0.8, "stroke_opacity": 0.3,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1},
        ).shift(LEFT * 3.2 + DOWN * 0.5)

        # 右侧: 不可对角化
        plane_r = NumberPlane(
            x_range=[-4, 4, 1], y_range=[-3, 3, 1],
            x_length=5, y_length=3.8,
            background_line_style={
                "stroke_color": GREEN_E, "stroke_width": 0.8, "stroke_opacity": 0.3,
            },
            axis_config={"stroke_color": GREY_B, "stroke_width": 1},
        ).shift(RIGHT * 3.2 + DOWN * 0.5)

        label_l = Text("Diagonalizable", font_size=22, color=GREEN_C).next_to(
            plane_l, UP, buff=0.15
        )
        label_r = Text("Non-diagonalizable", font_size=22, color=RED_B).next_to(
            plane_r, UP, buff=0.15
        )

        mat_l = MathTex(
            r"\begin{pmatrix}4&2\\1&3\end{pmatrix}",
            font_size=24,
        ).next_to(label_l, UP, buff=0.1)
        mat_r = MathTex(
            r"\begin{pmatrix}3&1\\0&3\end{pmatrix}",
            font_size=24,
        ).next_to(label_r, UP, buff=0.1)

        self.play(
            Create(plane_l), Create(plane_r),
            Write(label_l), Write(label_r),
            Write(mat_l), Write(mat_r),
            run_time=1.5,
        )

        # 网格变换对比
        identity = np.eye(2)
        grid_l_before = create_grid_lines(
            plane_l, identity, BLUE_D, num_lines=7,
            x_range=(-3, 3), y_range=(-2.5, 2.5), opacity=0.4,
        )
        grid_r_before = create_grid_lines(
            plane_r, identity, BLUE_D, num_lines=7,
            x_range=(-3, 3), y_range=(-2.5, 2.5), opacity=0.4,
        )

        grid_l_after = create_grid_lines(
            plane_l, A_good, GREEN_C, num_lines=7,
            x_range=(-3, 3), y_range=(-2.5, 2.5), opacity=0.5,
        )
        grid_r_after = create_grid_lines(
            plane_r, A_defective, RED_B, num_lines=7,
            x_range=(-3, 3), y_range=(-2.5, 2.5), opacity=0.5,
        )

        self.play(Create(grid_l_before), Create(grid_r_before), run_time=1)
        self.wait(0.5)

        self.play(
            Transform(grid_l_before, grid_l_after),
            Transform(grid_r_before, grid_r_after),
            run_time=2.5,
        )
        self.wait(0.5)

        # 特征方向标注
        # 左侧两个特征方向
        ev1_l = np.array([2, 1], dtype=float)
        ev2_l = np.array([-1, 1], dtype=float)
        ev1_l_n = ev1_l / np.linalg.norm(ev1_l)
        ev2_l_n = ev2_l / np.linalg.norm(ev2_l)

        dir_l1 = DashedLine(
            plane_l.c2p(*(- ev1_l_n * 2.5)),
            plane_l.c2p(*(ev1_l_n * 2.5)),
            color=YELLOW, stroke_width=2.5,
        )
        dir_l2 = DashedLine(
            plane_l.c2p(*(- ev2_l_n * 2.5)),
            plane_l.c2p(*(ev2_l_n * 2.5)),
            color=YELLOW, stroke_width=2.5,
        )
        dirs_l_label = Text("2 eigendirections", font_size=16,
                            color=YELLOW).next_to(plane_l, DOWN, buff=0.15)

        # 右侧一个特征方向
        dir_r1 = DashedLine(
            plane_r.c2p(-2.5, 0),
            plane_r.c2p(2.5, 0),
            color=YELLOW, stroke_width=2.5,
        )
        dirs_r_label = Text("Only 1 eigendirection", font_size=16,
                            color=RED_B).next_to(plane_r, DOWN, buff=0.15)

        self.play(
            Create(dir_l1), Create(dir_l2), Write(dirs_l_label),
            Create(dir_r1), Write(dirs_r_label),
            run_time=1.5,
        )
        self.wait(1)

        # 底部总结
        summary = Text(
            "Diagonalizable: enough eigenvectors to span the space",
            font_size=20, color=GREEN_A,
        ).to_edge(DOWN)
        self.play(FadeIn(summary), run_time=0.8)
        self.wait(2.5)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 4: 通过对角化计算矩阵幂次
# A^n = P D^n P^{-1}
# 用 A = [[4,2],[1,3]], D = [[5,0],[0,2]]
# ============================================================
class MatrixPowerViaDiagonalization(Scene):
    def construct(self):
        A = np.array([[4, 2], [1, 3]], dtype=float)
        P = np.array([[2, -1], [1, 1]], dtype=float)
        P_inv = np.linalg.inv(P)
        D = np.diag([5, 2])

        # ===== 阶段1: 引入公式 =====
        title = Text("Matrix Powers via Diagonalization", font_size=30).to_edge(UP)
        self.play(Write(title), run_time=1)

        # 核心公式推导
        eq1 = MathTex(r"A = PDP^{-1}", font_size=38)
        eq2 = MathTex(r"A^2 = (PDP^{-1})(PDP^{-1}) = PD^2P^{-1}", font_size=34)
        eq3 = MathTex(r"A^n = PD^nP^{-1}", font_size=42, color=YELLOW)

        eqs = VGroup(eq1, eq2, eq3).arrange(DOWN, buff=0.6).shift(UP * 0.3)

        self.play(Write(eq1), run_time=1.2)
        self.wait(0.5)
        self.play(Write(eq2), run_time=1.5)
        self.wait(0.5)
        self.play(Write(eq3), run_time=1.2)

        eq3_box = SurroundingRectangle(eq3, color=YELLOW, buff=0.15)
        self.play(Create(eq3_box), run_time=0.8)
        self.wait(1)

        # D^n 的计算
        dn_eq = MathTex(
            r"D^n = \begin{pmatrix} 5^n & 0 \\ 0 & 2^n \end{pmatrix}",
            font_size=36,
        ).next_to(eqs, DOWN, buff=0.6)
        dn_note = Text(
            "Diagonal powers: just raise each diagonal entry!",
            font_size=20, color=GREEN_A,
        ).next_to(dn_eq, DOWN, buff=0.2)

        self.play(Write(dn_eq), FadeIn(dn_note), run_time=1.5)
        self.wait(2)

        # 清理
        self.play(*[FadeOut(m) for m in self.mobjects if m != title], run_time=0.8)

        # ===== 阶段2: 几何展示 A, A^2, A^3 的效果 =====
        self.play(
            title.animate.become(
                Text("Geometric Effect of Matrix Powers", font_size=28).to_edge(UP)
            ),
            run_time=0.8,
        )

        # 三个子图
        planes = []
        grids_before = []
        grids_after = []
        power_labels = []

        x_positions = [-4, 0, 4]
        power_names = [r"A", r"A^2", r"A^3"]
        matrices = [A, np.linalg.matrix_power(A, 2), np.linalg.matrix_power(A, 3)]
        power_colors = [BLUE_C, GREEN_C, ORANGE]

        for idx in range(3):
            p = NumberPlane(
                x_range=[-3, 3, 1], y_range=[-2, 2, 1],
                x_length=3.5, y_length=2.8,
                background_line_style={
                    "stroke_color": GREY_E, "stroke_width": 0.6, "stroke_opacity": 0.25,
                },
                axis_config={"stroke_color": GREY_B, "stroke_width": 1},
            ).shift(RIGHT * x_positions[idx] + DOWN * 0.6)
            planes.append(p)

            lbl = MathTex(power_names[idx], font_size=32,
                          color=power_colors[idx]).next_to(p, UP, buff=0.15)
            power_labels.append(lbl)

            # 单位圆作为变换对象
            circle_pts_before = []
            circle_pts_after = []
            thetas = np.linspace(0, 2 * np.pi, 100)
            for th in thetas:
                pt = np.array([np.cos(th), np.sin(th)])
                circle_pts_before.append(p.c2p(*pt))
                transformed = matrices[idx] @ pt
                # 缩放以适应子图
                scale_factor = 1.0 / max(np.max(np.abs(matrices[idx])), 1)
                circle_pts_after.append(p.c2p(*(transformed * scale_factor)))

            circle_b = VMobject(color=BLUE_D, stroke_width=2, stroke_opacity=0.6)
            circle_b.set_points_smoothly(circle_pts_before)
            grids_before.append(circle_b)

            circle_a = VMobject(color=power_colors[idx], stroke_width=2.5)
            circle_a.set_points_smoothly(circle_pts_after)
            grids_after.append(circle_a)

        # 动画展示
        self.play(
            *[Create(p) for p in planes],
            *[Write(l) for l in power_labels],
            run_time=1.5,
        )
        self.play(
            *[Create(g) for g in grids_before],
            run_time=1,
        )
        self.wait(0.5)

        # 逐个变换
        for idx in range(3):
            self.play(
                Transform(grids_before[idx], grids_after[idx]),
                run_time=1.5,
            )
            self.wait(0.3)

        self.wait(1)

        # 在每个子图下方添加特征值幂次
        ev_pow_labels = [
            MathTex(r"5^1\!=\!5,\; 2^1\!=\!2", font_size=18, color=BLUE_C),
            MathTex(r"5^2\!=\!25,\; 2^2\!=\!4", font_size=18, color=GREEN_C),
            MathTex(r"5^3\!=\!125,\; 2^3\!=\!8", font_size=18, color=ORANGE),
        ]
        for idx in range(3):
            ev_pow_labels[idx].next_to(planes[idx], DOWN, buff=0.15)

        self.play(
            *[Write(l) for l in ev_pow_labels],
            run_time=1,
        )
        self.wait(1)

        note = Text(
            "Eigenvalues grow exponentially -> transformation stretches more",
            font_size=18, color=GREY_A,
        ).to_edge(DOWN)
        self.play(FadeIn(note), run_time=0.8)
        self.wait(1.5)

        # ===== 阶段3: 直接计算 vs 对角化计算的对比 =====
        self.play(*[FadeOut(m) for m in self.mobjects if m != title], run_time=0.8)
        self.play(
            title.animate.become(
                Text("Direct Computation vs Diagonalization", font_size=28).to_edge(UP)
            ),
            run_time=0.8,
        )

        # 左侧: 直接计算
        direct_title = Text("Direct: multiply n times", font_size=22,
                            color=RED_C).shift(LEFT * 3.2 + UP * 1.8)
        direct_lines = VGroup(
            MathTex(r"A^{10} = A \cdot A \cdot A \cdots A", font_size=26),
            MathTex(r"\text{10 matrix multiplications}", font_size=22, color=GREY_A),
            MathTex(r"\text{Each: } O(n^3) \text{ operations}", font_size=22,
                    color=GREY_A),
            MathTex(r"\text{Total: } O(n^3 \times k)", font_size=22, color=RED_C),
        ).arrange(DOWN, buff=0.3, aligned_edge=LEFT).next_to(direct_title, DOWN, buff=0.3)

        # 右侧: 对角化
        diag_title = Text("Diagonalization: compute once", font_size=22,
                          color=GREEN_C).shift(RIGHT * 3.2 + UP * 1.8)
        diag_lines = VGroup(
            MathTex(r"A^{10} = P D^{10} P^{-1}", font_size=26),
            MathTex(
                r"D^{10} = \begin{pmatrix} 5^{10} & 0 \\ 0 & 2^{10} \end{pmatrix}",
                font_size=24,
            ),
            MathTex(r"= \begin{pmatrix} 9765625 & 0 \\ 0 & 1024 \end{pmatrix}",
                    font_size=24, color=GREEN_C),
            MathTex(r"\text{Total: } O(n^3 + n \log k)", font_size=22,
                    color=GREEN_C),
        ).arrange(DOWN, buff=0.3, aligned_edge=LEFT).next_to(diag_title, DOWN, buff=0.3)

        # 分隔线
        divider = DashedLine(UP * 2, DOWN * 2.5, color=GREY,
                             stroke_width=1.5)

        self.play(
            Create(divider),
            Write(direct_title), Write(diag_title),
            run_time=1,
        )
        self.play(
            LaggedStart(*[Write(l) for l in direct_lines], lag_ratio=0.4),
            LaggedStart(*[Write(l) for l in diag_lines], lag_ratio=0.4),
            run_time=3,
        )
        self.wait(1.5)

        # 最终总结
        summary = VGroup(
            Text("Diagonalization: one-time cost, repeated benefit",
                 font_size=22, color=YELLOW),
            MathTex(
                r"A^n = P \begin{pmatrix} \lambda_1^n & & \\ & \ddots & \\"
                r" & & \lambda_m^n \end{pmatrix} P^{-1}",
                font_size=32,
            ),
        ).arrange(DOWN, buff=0.3).to_edge(DOWN).shift(UP * 0.2)

        summary_rect = SurroundingRectangle(summary, color=YELLOW, buff=0.2,
                                            fill_opacity=0.08)

        self.play(FadeIn(summary), Create(summary_rect), run_time=1.2)
        self.wait(3)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)
