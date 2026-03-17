"""
第8章：分块矩阵与块对角化 - Manim 可视化动画
===========================================

文件用途：为线性代数教学第8章制作3个可视化场景
输出格式：MP4 动画
运行命令：
  conda activate teaching_env && manim -ql --no_latex_cleanup ch8_block_matrix_visualizations.py BlockMatrixMultiplication
  conda activate teaching_env && manim -ql --no_latex_cleanup ch8_block_matrix_visualizations.py BlockDiagonalization
  conda activate teaching_env && manim -ql --no_latex_cleanup ch8_block_matrix_visualizations.py SchurComplement

作者：kyksj-1
"""

from manim import *
import numpy as np


# ============================================================
# 辅助函数
# ============================================================

def create_matrix_grid(entries, row_count, col_count, cell_size=0.6,
                       font_size=28, color=WHITE):
    """
    创建一个网格矩阵显示，返回 VGroup 包含：
    - grid_lines: 网格线
    - entries_tex: 每个单元格的文字
    - brackets: 左右方括号
    整体结构便于分块着色
    """
    group = VGroup()
    w = col_count * cell_size
    h = row_count * cell_size

    # 网格线（内部线）
    grid_lines = VGroup()
    for i in range(1, row_count):
        line = Line(
            LEFT * w / 2 + DOWN * (i * cell_size - h / 2),
            RIGHT * w / 2 + DOWN * (i * cell_size - h / 2),
            stroke_width=1, color=GREY_B, stroke_opacity=0.5,
        )
        grid_lines.add(line)
    for j in range(1, col_count):
        line = Line(
            UP * h / 2 + RIGHT * (j * cell_size - w / 2),
            DOWN * h / 2 + RIGHT * (j * cell_size - w / 2),
            stroke_width=1, color=GREY_B, stroke_opacity=0.5,
        )
        grid_lines.add(line)
    group.add(grid_lines)

    # 单元格文字
    entries_tex = VGroup()
    for i in range(row_count):
        for j in range(col_count):
            idx = i * col_count + j
            text = entries[idx] if idx < len(entries) else ""
            t = MathTex(str(text), font_size=font_size, color=color)
            # 定位：左上角为原点
            x_pos = (j + 0.5) * cell_size - w / 2
            y_pos = h / 2 - (i + 0.5) * cell_size
            t.move_to(np.array([x_pos, y_pos, 0]))
            entries_tex.add(t)
    group.add(entries_tex)

    # 方括号
    left_bracket = MathTex(r"\Bigg[", font_size=int(font_size * row_count * 1.1))
    left_bracket.next_to(group, LEFT, buff=0.05)
    right_bracket = MathTex(r"\Bigg]", font_size=int(font_size * row_count * 1.1))
    right_bracket.next_to(group, RIGHT, buff=0.05)
    brackets = VGroup(left_bracket, right_bracket)
    group.add(brackets)

    return group, grid_lines, entries_tex, brackets


def create_block_highlight(row_start, row_end, col_start, col_end,
                           row_count, col_count, cell_size=0.6,
                           color=BLUE, opacity=0.25, stroke_width=2.5):
    """
    为矩阵网格的指定分块区域创建高亮矩形
    行列索引从 0 开始
    """
    w = col_count * cell_size
    h = row_count * cell_size

    x_left = col_start * cell_size - w / 2
    x_right = (col_end + 1) * cell_size - w / 2
    y_top = h / 2 - row_start * cell_size
    y_bottom = h / 2 - (row_end + 1) * cell_size

    rect = Rectangle(
        width=x_right - x_left,
        height=y_top - y_bottom,
        fill_color=color,
        fill_opacity=opacity,
        stroke_color=color,
        stroke_width=stroke_width,
        stroke_opacity=0.9,
    )
    rect.move_to(np.array([(x_left + x_right) / 2,
                            (y_top + y_bottom) / 2, 0]))
    return rect


# ============================================================
# Scene 1: 分块矩阵乘法
# 展示 2x2 分块结构的矩阵乘法，"块乘块"的过程
# M (4x4) = A (4x4) * B (4x4)，各自分成 2x2 的分块
# ============================================================
class BlockMatrixMultiplication(Scene):
    def construct(self):
        # ---------- 阶段1: 标题引入 ----------
        title = Text("Block Matrix Multiplication", font_size=32).to_edge(UP)
        self.play(Write(title), run_time=1)

        subtitle = Text(
            "Treat blocks as single elements",
            font_size=22, color=GREY_A,
        ).next_to(title, DOWN, buff=0.25)
        self.play(FadeIn(subtitle), run_time=0.8)
        self.wait(1)

        # ---------- 阶段2: 展示分块矩阵的公式 ----------
        # 分块矩阵乘法公式
        formula = MathTex(
            r"\begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}"
            r"\begin{pmatrix} B_{11} & B_{12} \\ B_{21} & B_{22} \end{pmatrix}"
            r"="
            r"\begin{pmatrix} C_{11} & C_{12} \\ C_{21} & C_{22} \end{pmatrix}",
            font_size=30,
        ).shift(UP * 0.8)

        self.play(Write(formula), run_time=2)
        self.wait(1)

        # 逐项公式
        detail = MathTex(
            r"C_{ij} = \sum_k A_{ik} B_{kj}",
            font_size=32, color=YELLOW,
        ).next_to(formula, DOWN, buff=0.5)
        detail_box = SurroundingRectangle(detail, color=YELLOW, buff=0.12,
                                          fill_opacity=0.08)
        self.play(Write(detail), Create(detail_box), run_time=1.2)
        self.wait(1.5)

        # 清理公式，准备进入具体数值演示
        self.play(
            FadeOut(formula), FadeOut(detail), FadeOut(detail_box),
            FadeOut(subtitle),
            run_time=0.8,
        )

        # ---------- 阶段3: 构建具体的 4x4 分块矩阵 ----------
        # 矩阵 A (4x4)
        # 分块: A11=[[1,2],[3,4]], A12=[[0,1],[1,0]], A21=[[1,0],[0,1]], A22=[[2,1],[1,2]]
        a_entries = [
            1, 2, 0, 1,
            3, 4, 1, 0,
            1, 0, 2, 1,
            0, 1, 1, 2,
        ]
        # 矩阵 B (4x4)
        b_entries = [
            1, 0, 1, 1,
            0, 1, 0, 1,
            2, 1, 1, 0,
            1, 0, 0, 1,
        ]

        cell_size = 0.55

        # --- 矩阵 A ---
        a_group, a_grid, a_tex, a_brackets = create_matrix_grid(
            a_entries, 4, 4, cell_size=cell_size, font_size=26,
        )
        a_group.shift(LEFT * 3.8 + DOWN * 0.3)
        a_label = MathTex("A", font_size=32, color=BLUE_C).next_to(a_group, UP, buff=0.2)

        # --- 矩阵 B ---
        b_group, b_grid, b_tex, b_brackets = create_matrix_grid(
            b_entries, 4, 4, cell_size=cell_size, font_size=26,
        )
        b_group.shift(LEFT * 0.2 + DOWN * 0.3)
        b_label = MathTex("B", font_size=32, color=GREEN_C).next_to(b_group, UP, buff=0.2)

        # 乘号
        times_sign = MathTex(r"\times", font_size=36).move_to(
            (a_group.get_right() + b_group.get_left()) / 2
        )

        # 等号
        eq_sign = MathTex(r"=", font_size=36).next_to(b_group, RIGHT, buff=0.35)

        # --- 结果矩阵 C (先用占位符) ---
        # 计算实际结果
        A_mat = np.array(a_entries).reshape(4, 4)
        B_mat = np.array(b_entries).reshape(4, 4)
        C_mat = A_mat @ B_mat
        c_entries = C_mat.flatten().tolist()

        # 先创建空白结果矩阵
        c_entries_q = ["?"] * 16
        c_group, c_grid, c_tex, c_brackets = create_matrix_grid(
            c_entries_q, 4, 4, cell_size=cell_size, font_size=26, color=GREY_B,
        )
        c_group.shift(RIGHT * 4.0 + DOWN * 0.3)
        c_label = MathTex("C = AB", font_size=28, color=RED_B).next_to(c_group, UP, buff=0.2)

        # 播放矩阵出现
        self.play(
            FadeIn(a_group), Write(a_label),
            FadeIn(b_group), Write(b_label),
            Write(times_sign), Write(eq_sign),
            FadeIn(c_group), Write(c_label),
            run_time=1.5,
        )
        self.wait(0.5)

        # ---------- 阶段4: 用颜色标示分块 ----------
        # 定义分块颜色
        block_colors = {
            (0, 0): BLUE,       # 左上
            (0, 1): TEAL,       # 右上
            (1, 0): ORANGE,     # 左下
            (1, 1): RED_B,      # 右下
        }

        # 为矩阵A着色分块
        a_highlights = VGroup()
        for (bi, bj), col in block_colors.items():
            r_start, r_end = bi * 2, bi * 2 + 1
            c_start, c_end = bj * 2, bj * 2 + 1
            rect = create_block_highlight(
                r_start, r_end, c_start, c_end, 4, 4,
                cell_size=cell_size, color=col, opacity=0.2,
            )
            rect.shift(a_group.get_center() - np.array([0, 0, 0]))
            # 修正位置：grid_lines的中心就是a_group的前两个子元素的中心
            rect.move_to(a_grid.get_center() + rect.get_center())
            a_highlights.add(rect)

        # 为矩阵B着色分块
        b_highlights = VGroup()
        for (bi, bj), col in block_colors.items():
            r_start, r_end = bi * 2, bi * 2 + 1
            c_start, c_end = bj * 2, bj * 2 + 1
            rect = create_block_highlight(
                r_start, r_end, c_start, c_end, 4, 4,
                cell_size=cell_size, color=col, opacity=0.2,
            )
            rect.move_to(b_grid.get_center() + rect.get_center())
            b_highlights.add(rect)

        block_text = Text(
            "Color-coded 2x2 blocks",
            font_size=20, color=GREY_A,
        ).to_edge(DOWN)

        self.play(
            LaggedStart(*[FadeIn(h) for h in a_highlights], lag_ratio=0.1),
            LaggedStart(*[FadeIn(h) for h in b_highlights], lag_ratio=0.1),
            FadeIn(block_text),
            run_time=1.5,
        )
        self.wait(1.5)

        # ---------- 阶段5: 演示 C11 = A11*B11 + A12*B21 ----------
        self.play(FadeOut(block_text), run_time=0.3)

        # 计算说明
        c11_formula = MathTex(
            r"C_{11} = A_{11} B_{11} + A_{12} B_{21}",
            font_size=28, color=YELLOW,
        ).to_edge(DOWN).shift(UP * 0.3)

        self.play(Write(c11_formula), run_time=1)
        self.wait(0.5)

        # 高亮 A11 (蓝) 和 B11 (蓝) 的贡献
        # 先让其他分块变暗
        # 高亮 A 的第一行分块
        a_row_highlight = VGroup()
        for bj in range(2):
            col = block_colors[(0, bj)]
            r_start, r_end = 0, 1
            c_start, c_end = bj * 2, bj * 2 + 1
            rect = create_block_highlight(
                r_start, r_end, c_start, c_end, 4, 4,
                cell_size=cell_size, color=col, opacity=0.35, stroke_width=3.5,
            )
            rect.move_to(a_grid.get_center() + rect.get_center())
            a_row_highlight.add(rect)

        # 高亮 B 的第一列分块
        b_col_highlight = VGroup()
        for bi in range(2):
            col = block_colors[(bi, 0)]
            r_start, r_end = bi * 2, bi * 2 + 1
            c_start, c_end = 0, 1
            rect = create_block_highlight(
                r_start, r_end, c_start, c_end, 4, 4,
                cell_size=cell_size, color=col, opacity=0.35, stroke_width=3.5,
            )
            rect.move_to(b_grid.get_center() + rect.get_center())
            b_col_highlight.add(rect)

        # 高亮结果 C11 区域
        c11_highlight = create_block_highlight(
            0, 1, 0, 1, 4, 4,
            cell_size=cell_size, color=YELLOW, opacity=0.3, stroke_width=3.5,
        )
        c11_highlight.move_to(c_grid.get_center() + c11_highlight.get_center())

        self.play(
            FadeOut(a_highlights), FadeOut(b_highlights),
            run_time=0.5,
        )
        self.play(
            LaggedStart(*[FadeIn(h) for h in a_row_highlight], lag_ratio=0.15),
            LaggedStart(*[FadeIn(h) for h in b_col_highlight], lag_ratio=0.15),
            FadeIn(c11_highlight),
            run_time=1.2,
        )
        self.wait(1)

        # 第一项: A11 * B11
        step1_text = MathTex(
            r"A_{11} B_{11} = "
            r"\begin{pmatrix}1&2\\3&4\end{pmatrix}"
            r"\begin{pmatrix}1&0\\0&1\end{pmatrix}"
            r"= \begin{pmatrix}1&2\\3&4\end{pmatrix}",
            font_size=24,
        ).next_to(c11_formula, DOWN, buff=0.3)

        self.play(Write(step1_text), run_time=1.5)
        self.wait(1)

        # 第二项: A12 * B21
        step2_text = MathTex(
            r"A_{12} B_{21} = "
            r"\begin{pmatrix}0&1\\1&0\end{pmatrix}"
            r"\begin{pmatrix}2&1\\1&0\end{pmatrix}"
            r"= \begin{pmatrix}1&0\\2&1\end{pmatrix}",
            font_size=24,
        ).next_to(step1_text, DOWN, buff=0.2)

        self.play(Write(step2_text), run_time=1.5)
        self.wait(1)

        # 清理中间步骤
        self.play(FadeOut(step1_text), FadeOut(step2_text), run_time=0.5)

        # 显示 C11 的实际结果
        result_text = MathTex(
            r"C_{11} = \begin{pmatrix}1&2\\3&4\end{pmatrix} + \begin{pmatrix}1&0\\2&1\end{pmatrix}"
            r"= \begin{pmatrix}2&2\\5&5\end{pmatrix}",
            font_size=26, color=YELLOW,
        ).to_edge(DOWN).shift(UP * 0.3)

        self.play(
            ReplacementTransform(c11_formula, result_text),
            run_time=1.2,
        )
        self.wait(1)

        # 填入 C11 的数值
        c11_values = [int(C_mat[i][j]) for i in range(2) for j in range(2)]
        new_c11_entries = VGroup()
        w = 4 * cell_size
        h = 4 * cell_size
        for i in range(2):
            for j in range(2):
                val = c11_values[i * 2 + j]
                t = MathTex(str(val), font_size=26, color=YELLOW)
                x_pos = (j + 0.5) * cell_size - w / 2
                y_pos = h / 2 - (i + 0.5) * cell_size
                t.move_to(c_grid.get_center() + np.array([x_pos, y_pos, 0]))
                new_c11_entries.add(t)

        # 替换 C 矩阵中 C11 位置的 "?"
        old_q = VGroup(*[c_tex[i * 4 + j] for i in range(2) for j in range(2)])
        self.play(
            ReplacementTransform(old_q, new_c11_entries),
            run_time=1,
        )
        self.wait(1)

        # ---------- 阶段6: 快速填入其余分块 ----------
        self.play(
            FadeOut(a_row_highlight), FadeOut(b_col_highlight),
            FadeOut(c11_highlight), FadeOut(result_text),
            run_time=0.5,
        )

        remaining_text = Text(
            "Similarly compute C12, C21, C22 ...",
            font_size=22, color=GREY_A,
        ).to_edge(DOWN)
        self.play(FadeIn(remaining_text), run_time=0.5)
        self.wait(1)

        # 填入所有剩余数值
        all_new_entries = VGroup()
        w = 4 * cell_size
        h = 4 * cell_size
        for i in range(4):
            for j in range(4):
                if i < 2 and j < 2:
                    continue  # C11 已填入
                val = int(C_mat[i][j])
                t = MathTex(str(val), font_size=26, color=WHITE)
                x_pos = (j + 0.5) * cell_size - w / 2
                y_pos = h / 2 - (i + 0.5) * cell_size
                t.move_to(c_grid.get_center() + np.array([x_pos, y_pos, 0]))
                all_new_entries.add(t)

        remaining_old = VGroup(*[
            c_tex[i * 4 + j]
            for i in range(4) for j in range(4)
            if not (i < 2 and j < 2)
        ])

        # 用分块颜色闪烁显示各个块
        c12_hl = create_block_highlight(0, 1, 2, 3, 4, 4,
                                        cell_size=cell_size, color=TEAL, opacity=0.3)
        c12_hl.move_to(c_grid.get_center() + c12_hl.get_center())
        c21_hl = create_block_highlight(2, 3, 0, 1, 4, 4,
                                        cell_size=cell_size, color=ORANGE, opacity=0.3)
        c21_hl.move_to(c_grid.get_center() + c21_hl.get_center())
        c22_hl = create_block_highlight(2, 3, 2, 3, 4, 4,
                                        cell_size=cell_size, color=RED_B, opacity=0.3)
        c22_hl.move_to(c_grid.get_center() + c22_hl.get_center())

        self.play(
            FadeIn(c12_hl), FadeIn(c21_hl), FadeIn(c22_hl),
            ReplacementTransform(remaining_old, all_new_entries),
            run_time=1.5,
        )
        self.wait(1)

        # ---------- 阶段7: 总结 ----------
        self.play(FadeOut(remaining_text), run_time=0.3)
        self.play(FadeOut(c12_hl), FadeOut(c21_hl), FadeOut(c22_hl), run_time=0.3)

        summary = Text(
            "Block multiplication mirrors scalar multiplication: row-of-blocks times column-of-blocks",
            font_size=20, color=GREEN_A,
        ).to_edge(DOWN)
        self.play(FadeIn(summary), run_time=1)
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 2: 块对角化——"分而治之"
# 展示一个 4x4 矩阵通过相似变换变为块对角形式
# 每个块独立处理，标注各块的特征值
# ============================================================
class BlockDiagonalization(Scene):
    def construct(self):
        # ---------- 阶段1: 标题和核心思想 ----------
        title = Text("Block Diagonalization: Divide and Conquer",
                      font_size=30).to_edge(UP)
        self.play(Write(title), run_time=1)

        # 核心等式
        main_eq = MathTex(
            r"P^{-1} A P = \begin{pmatrix} B_1 & 0 \\ 0 & B_2 \end{pmatrix}",
            font_size=36, color=YELLOW,
        ).shift(UP * 1.5)
        main_eq_box = SurroundingRectangle(main_eq, color=YELLOW, buff=0.15,
                                            fill_opacity=0.06)

        self.play(Write(main_eq), Create(main_eq_box), run_time=1.5)
        self.wait(1)

        idea_text = Text(
            "Each block can be analyzed independently!",
            font_size=24, color=TEAL_A,
        ).next_to(main_eq_box, DOWN, buff=0.4)
        self.play(FadeIn(idea_text), run_time=0.8)
        self.wait(1.5)

        self.play(
            FadeOut(main_eq), FadeOut(main_eq_box), FadeOut(idea_text),
            run_time=0.8,
        )

        # ---------- 阶段2: 展示原始矩阵 ----------
        self.play(
            title.animate.become(
                Text("Example: 4x4 Block Diagonalization", font_size=28).to_edge(UP)
            ),
            run_time=0.6,
        )

        # 原始 4x4 矩阵（有结构但不明显）
        # A = P @ diag(B1, B2) @ P^{-1}
        # B1 = [[3,1],[0,2]], B2 = [[1,-1],[1,1]]
        original_mat = MathTex(
            r"A = \begin{pmatrix}"
            r"3 & 1 & 0 & 0 \\"
            r"0 & 2 & 0 & 0 \\"
            r"0 & 0 & 1 & -1 \\"
            r"0 & 0 & 1 & 1"
            r"\end{pmatrix}",
            font_size=32,
        ).shift(LEFT * 3 + DOWN * 0.3)
        orig_box = SurroundingRectangle(original_mat, color=BLUE, buff=0.15,
                                         fill_opacity=0.05, stroke_width=1.5)

        self.play(Write(original_mat), Create(orig_box), run_time=1.5)
        self.wait(1)

        # 标注"观察矩阵有两个独立的2x2块"
        observe_text = Text(
            "Observe: A is already block diagonal!",
            font_size=22, color=GREY_A,
        ).to_edge(DOWN)
        self.play(FadeIn(observe_text), run_time=0.6)
        self.wait(1)

        # ---------- 阶段3: 高亮两个独立分块 ----------
        # 高亮左上块 B1
        b1_rect = Rectangle(
            width=2.0, height=1.2,
            stroke_color=BLUE_C, stroke_width=3,
            fill_color=BLUE_C, fill_opacity=0.15,
        )
        # 手动定位到 A 矩阵的左上 2x2 子块
        b1_rect.move_to(original_mat.get_center() + UP * 0.52 + LEFT * 0.55)

        # 高亮右下块 B2
        b2_rect = Rectangle(
            width=2.0, height=1.2,
            stroke_color=RED_B, stroke_width=3,
            fill_color=RED_B, fill_opacity=0.15,
        )
        b2_rect.move_to(original_mat.get_center() + DOWN * 0.55 + RIGHT * 0.55)

        b1_label = MathTex(r"B_1", font_size=28, color=BLUE_C).next_to(b1_rect, UP, buff=0.1)
        b2_label = MathTex(r"B_2", font_size=28, color=RED_B).next_to(b2_rect, DOWN, buff=0.1)

        self.play(
            Create(b1_rect), Write(b1_label),
            Create(b2_rect), Write(b2_label),
            run_time=1.2,
        )
        self.wait(1.5)

        # ---------- 阶段4: "分而治之"——分别显示每个块 ----------
        self.play(FadeOut(observe_text), run_time=0.3)

        # 箭头指向右边
        arrow = Arrow(
            original_mat.get_right() + RIGHT * 0.3,
            original_mat.get_right() + RIGHT * 1.8,
            buff=0, color=WHITE, stroke_width=3,
        )
        arrow_label = Text("Divide", font_size=20, color=WHITE).next_to(arrow, UP, buff=0.1)

        self.play(GrowArrow(arrow), FadeIn(arrow_label), run_time=0.8)

        # 右侧展示两个独立子块
        b1_display = MathTex(
            r"B_1 = \begin{pmatrix} 3 & 1 \\ 0 & 2 \end{pmatrix}",
            font_size=30, color=BLUE_C,
        ).shift(RIGHT * 3 + UP * 1.2)
        b1_box = SurroundingRectangle(b1_display, color=BLUE_C, buff=0.12,
                                       fill_opacity=0.06)

        b2_display = MathTex(
            r"B_2 = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix}",
            font_size=30, color=RED_B,
        ).shift(RIGHT * 3 + DOWN * 1.2)
        b2_box = SurroundingRectangle(b2_display, color=RED_B, buff=0.12,
                                       fill_opacity=0.06)

        self.play(
            Write(b1_display), Create(b1_box),
            Write(b2_display), Create(b2_box),
            run_time=1.5,
        )
        self.wait(1)

        # ---------- 阶段5: 标注每个块的特征值 ----------
        # B1 = [[3,1],[0,2]]: 上三角，特征值 = 3, 2
        b1_eigenvals = MathTex(
            r"\lambda(B_1) = \{3, 2\}",
            font_size=24, color=YELLOW,
        ).next_to(b1_box, DOWN, buff=0.2)

        # B2 = [[1,-1],[1,1]]: 特征多项式 (1-lam)^2+1=0 => lam = 1+/-i
        b2_eigenvals = MathTex(
            r"\lambda(B_2) = \{1+i, 1-i\}",
            font_size=24, color=YELLOW,
        ).next_to(b2_box, DOWN, buff=0.2)

        self.play(Write(b1_eigenvals), Write(b2_eigenvals), run_time=1.2)
        self.wait(1)

        # 全局特征值 = 两块的并集
        global_eig = MathTex(
            r"\lambda(A) = \lambda(B_1) \cup \lambda(B_2) = \{3, 2, 1+i, 1-i\}",
            font_size=28, color=GREEN_A,
        ).to_edge(DOWN).shift(UP * 0.2)
        global_box = SurroundingRectangle(global_eig, color=GREEN_A, buff=0.1,
                                          fill_opacity=0.06)

        self.play(Write(global_eig), Create(global_box), run_time=1.2)
        self.wait(1.5)

        # ---------- 阶段6: 强调"分而治之"效益 ----------
        self.play(
            FadeOut(b1_eigenvals), FadeOut(b2_eigenvals),
            FadeOut(global_eig), FadeOut(global_box),
            run_time=0.5,
        )

        # 复杂度对比
        complexity_title = Text("Computational Benefit", font_size=24,
                                color=TEAL_A).to_edge(DOWN).shift(UP * 1.5)

        cost_full = MathTex(
            r"\text{Full } 4 \times 4: \quad O(4^3) = O(64)",
            font_size=26,
        ).next_to(complexity_title, DOWN, buff=0.3)

        cost_block = MathTex(
            r"\text{Two } 2 \times 2: \quad 2 \times O(2^3) = O(16)",
            font_size=26, color=GREEN_A,
        ).next_to(cost_full, DOWN, buff=0.2)

        speedup = MathTex(
            r"\text{Speedup: } 4\times",
            font_size=28, color=YELLOW,
        ).next_to(cost_block, DOWN, buff=0.3)

        self.play(
            Write(complexity_title),
            run_time=0.8,
        )
        self.play(Write(cost_full), run_time=1)
        self.play(Write(cost_block), run_time=1)
        self.play(Write(speedup), run_time=0.8)
        self.wait(2)

        # ---------- 阶段7: 总结 ----------
        self.play(
            FadeOut(complexity_title), FadeOut(cost_full),
            FadeOut(cost_block), FadeOut(speedup),
            run_time=0.5,
        )

        summary = VGroup(
            Text("Block diagonal structure:", font_size=24, color=WHITE),
            Text("1. Eigenvalues = union of block eigenvalues", font_size=22, color=GREY_A),
            Text("2. Determinant = product of block determinants", font_size=22, color=GREY_A),
            Text("3. Inverse = block diagonal of block inverses", font_size=22, color=GREY_A),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT).to_edge(DOWN).shift(UP * 0.2)

        self.play(
            LaggedStart(*[FadeIn(s) for s in summary], lag_ratio=0.3),
            run_time=2,
        )
        self.wait(2.5)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)


# ============================================================
# Scene 3: Schur 补与分块消元
# 2x2 分块矩阵 M = [[A, B],[C, D]]
# 通过分块行变换消去 C，得到 Schur 补 S = D - C A^{-1} B
# ============================================================
class SchurComplement(Scene):
    def construct(self):
        # ---------- 阶段1: 标题和问题引入 ----------
        title = Text("Schur Complement and Block Elimination",
                      font_size=30).to_edge(UP)
        self.play(Write(title), run_time=1)

        # 问题：给定分块矩阵，如何做分块"高斯消元"？
        problem_text = Text(
            "Given a 2x2 block matrix, how to perform block Gaussian elimination?",
            font_size=22, color=GREY_A,
        ).next_to(title, DOWN, buff=0.3)
        self.play(FadeIn(problem_text), run_time=0.8)
        self.wait(1)

        # ---------- 阶段2: 展示分块矩阵 M ----------
        self.play(FadeOut(problem_text), run_time=0.3)

        # 原始分块矩阵
        block_M = MathTex(
            r"M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}",
            font_size=38,
        ).shift(UP * 1.2)

        self.play(Write(block_M), run_time=1.2)
        self.wait(1)

        # ---------- 阶段3: 分块消元思路——类比标量消元 ----------
        # 标量情况回顾
        scalar_title = Text("Recall: Scalar Gaussian Elimination",
                            font_size=22, color=TEAL_A).shift(DOWN * 0.3)

        scalar_eq = MathTex(
            r"\begin{pmatrix} a & b \\ c & d \end{pmatrix}"
            r"\overset{R_2 - \frac{c}{a} R_1}{\longrightarrow}"
            r"\begin{pmatrix} a & b \\ 0 & d - \frac{c}{a} b \end{pmatrix}",
            font_size=30,
        ).next_to(scalar_title, DOWN, buff=0.4)

        self.play(Write(scalar_title), run_time=0.8)
        self.play(Write(scalar_eq), run_time=1.5)
        self.wait(1.5)

        # ---------- 阶段4: 推广到分块——核心公式 ----------
        self.play(FadeOut(scalar_title), FadeOut(scalar_eq), run_time=0.5)

        # 分块消元：左乘 [[I, 0],[-CA^{-1}, I]]
        elim_matrix = MathTex(
            r"\begin{pmatrix} I & 0 \\ -CA^{-1} & I \end{pmatrix}",
            font_size=32, color=TEAL_B,
        ).shift(LEFT * 3.5 + DOWN * 0.5)

        elim_times = MathTex(r"\cdot", font_size=32).next_to(elim_matrix, RIGHT, buff=0.15)

        original_M = MathTex(
            r"\begin{pmatrix} A & B \\ C & D \end{pmatrix}",
            font_size=32,
        ).next_to(elim_times, RIGHT, buff=0.15)

        eq_sign = MathTex(r"=", font_size=32).next_to(original_M, RIGHT, buff=0.2)

        result_M = MathTex(
            r"\begin{pmatrix} A & B \\ 0 & D - CA^{-1}B \end{pmatrix}",
            font_size=32, color=YELLOW,
        ).next_to(eq_sign, RIGHT, buff=0.2)

        step_label = Text("Block row operation: R2 <- R2 - C A^(-1) R1",
                          font_size=20, color=GREY_A).to_edge(DOWN)

        self.play(
            Write(elim_matrix),
            Write(elim_times),
            Write(original_M),
            run_time=1.5,
        )
        self.play(
            Write(eq_sign),
            Write(result_M),
            FadeIn(step_label),
            run_time=1.5,
        )
        self.wait(2)

        # ---------- 阶段5: 突出 Schur 补 ----------
        self.play(FadeOut(step_label), run_time=0.3)

        # Schur 补定义
        schur_def = MathTex(
            r"S = D - C A^{-1} B",
            font_size=38, color=RED_B,
        ).shift(DOWN * 2.2)
        schur_box = SurroundingRectangle(schur_def, color=RED_B, buff=0.15,
                                         fill_opacity=0.08)
        schur_label = Text("Schur Complement of A in M",
                           font_size=22, color=RED_B).next_to(schur_box, DOWN, buff=0.2)

        self.play(
            Write(schur_def), Create(schur_box),
            Write(schur_label),
            run_time=1.5,
        )
        self.wait(2)

        # ---------- 阶段6: 清理并展示具体数值例子 ----------
        self.play(
            FadeOut(block_M), FadeOut(elim_matrix), FadeOut(elim_times),
            FadeOut(original_M), FadeOut(eq_sign), FadeOut(result_M),
            FadeOut(schur_def), FadeOut(schur_box), FadeOut(schur_label),
            run_time=0.8,
        )

        self.play(
            title.animate.become(
                Text("Numerical Example", font_size=28).to_edge(UP)
            ),
            run_time=0.6,
        )

        # 具体矩阵
        # M = [[4,2,1,0],[1,3,0,1],[2,0,5,1],[0,1,1,3]]
        # A = [[4,2],[1,3]], B = [[1,0],[0,1]], C = [[2,0],[0,1]], D = [[5,1],[1,3]]
        num_M = MathTex(
            r"M = \begin{pmatrix}"
            r"4 & 2 & 1 & 0 \\"
            r"1 & 3 & 0 & 1 \\"
            r"2 & 0 & 5 & 1 \\"
            r"0 & 1 & 1 & 3"
            r"\end{pmatrix}",
            font_size=30,
        ).shift(UP * 1.5)

        self.play(Write(num_M), run_time=1.5)
        self.wait(0.5)

        # 标注分块
        blocks_text = MathTex(
            r"A = \begin{pmatrix}4&2\\1&3\end{pmatrix},\ "
            r"B = \begin{pmatrix}1&0\\0&1\end{pmatrix},\ "
            r"C = \begin{pmatrix}2&0\\0&1\end{pmatrix},\ "
            r"D = \begin{pmatrix}5&1\\1&3\end{pmatrix}",
            font_size=24,
        ).next_to(num_M, DOWN, buff=0.5)

        self.play(Write(blocks_text), run_time=1.5)
        self.wait(1)

        # ---------- 阶段7: 计算 Schur 补 ----------
        # A^{-1} = 1/10 * [[3,-2],[-1,4]]
        step_Ainv = MathTex(
            r"A^{-1} = \frac{1}{10}\begin{pmatrix}3&-2\\-1&4\end{pmatrix}",
            font_size=26,
        ).shift(DOWN * 0.6)

        self.play(Write(step_Ainv), run_time=1.2)
        self.wait(0.5)

        # C A^{-1} B
        step_CAinvB = MathTex(
            r"CA^{-1}B = \begin{pmatrix}2&0\\0&1\end{pmatrix}"
            r"\frac{1}{10}\begin{pmatrix}3&-2\\-1&4\end{pmatrix}"
            r"\begin{pmatrix}1&0\\0&1\end{pmatrix}"
            r"= \frac{1}{10}\begin{pmatrix}6&-4\\-1&4\end{pmatrix}",
            font_size=22,
        ).next_to(step_Ainv, DOWN, buff=0.35)

        self.play(Write(step_CAinvB), run_time=1.5)
        self.wait(0.5)

        # S = D - C A^{-1} B
        # S = [[5,1],[1,3]] - 0.1*[[6,-4],[-1,4]]
        #   = [[5-0.6, 1+0.4],[1+0.1, 3-0.4]]
        #   = [[4.4, 1.4],[1.1, 2.6]]
        step_S = MathTex(
            r"S = D - CA^{-1}B = \begin{pmatrix}5&1\\1&3\end{pmatrix}"
            r"- \frac{1}{10}\begin{pmatrix}6&-4\\-1&4\end{pmatrix}"
            r"= \begin{pmatrix}4.4&1.4\\1.1&2.6\end{pmatrix}",
            font_size=24, color=RED_B,
        ).next_to(step_CAinvB, DOWN, buff=0.35)
        step_S_box = SurroundingRectangle(step_S, color=RED_B, buff=0.1,
                                          fill_opacity=0.06)

        self.play(Write(step_S), Create(step_S_box), run_time=1.5)
        self.wait(1.5)

        # ---------- 阶段8: Schur 补的应用意义 ----------
        self.play(
            FadeOut(num_M), FadeOut(blocks_text),
            FadeOut(step_Ainv), FadeOut(step_CAinvB),
            FadeOut(step_S), FadeOut(step_S_box),
            run_time=0.8,
        )

        self.play(
            title.animate.become(
                Text("Why Schur Complement Matters", font_size=28).to_edge(UP)
            ),
            run_time=0.6,
        )

        # 关键性质列表
        prop_colors = [YELLOW, BLUE_C, GREEN_C]
        properties = VGroup(
            MathTex(
                r"\det(M) = \det(A) \cdot \det(S)",
                font_size=30, color=YELLOW,
            ),
            MathTex(
                r"M \succ 0 \Longleftrightarrow A \succ 0 \text{ and } S \succ 0",
                font_size=28, color=BLUE_C,
            ),
            Text(
                "Solving Mx = b reduces to two smaller systems",
                font_size=24, color=GREEN_C,
            ),
        ).arrange(DOWN, buff=0.6).shift(UP * 0.3)

        for i, prop in enumerate(properties):
            box = SurroundingRectangle(prop, color=prop_colors[i], buff=0.12,
                                       fill_opacity=0.05, stroke_width=1.5)
            self.play(Write(prop), Create(box), run_time=1.2)
            self.wait(1)

        self.wait(1)

        # ---------- 阶段9: 行列式验证 ----------
        det_verify = MathTex(
            r"\det(M) = \det(A) \cdot \det(S) = 10 \times (4.4 \times 2.6 - 1.4 \times 1.1)"
            r"= 10 \times 9.9 = 99",
            font_size=24, color=GREY_A,
        ).to_edge(DOWN).shift(UP * 0.2)

        self.play(FadeIn(det_verify), run_time=1)
        self.wait(2)

        # ---------- 阶段10: 总结 ----------
        summary = Text(
            "Schur complement: reduces block systems to smaller problems",
            font_size=22, color=GREEN_A,
        ).to_edge(DOWN)

        self.play(
            FadeOut(det_verify),
            FadeIn(summary),
            run_time=1,
        )
        self.wait(2)
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=1.5)
