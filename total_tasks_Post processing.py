# Author  : xuzhou
# Date    : 2019-08-15
"""
数据处理, 包含了所有功能
"""
import os
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm


class DrawSmoothCurves:
    """绘制平滑曲线类。"""

    @staticmethod
    def bezier_curve(p0, p1, p2, p3, inserted):
        """三阶贝塞尔曲线。

        :param p0: 点坐标, tuple、list或numpy.ndarray类型
        :param p1: 点坐标, tuple、list或numpy.ndarray类型
        :param p2: 点坐标, tuple、list或numpy.ndarray类型
        :param p3: 点坐标, tuple、list或numpy.ndarray类型
        :param inserted: p0和p3之间插值的数量
        :return: 插值完成后的点
        """
        # 对 p0, p1, p2, p3 进行类型检查
        if not isinstance(p0, (tuple, list, np.ndarray)):
            raise ValueError("p0点坐标应该是元组、列表或numpy数组类型")
        if not isinstance(p1, (tuple, list, np.ndarray)):
            raise ValueError("p1点坐标应该是元组、列表或numpy数组类型")
        if not isinstance(p2, (tuple, list, np.ndarray)):
            raise ValueError("p2点坐标应该是元组、列表或numpy数组类型")
        if not isinstance(p3, (tuple, list, np.ndarray)):
            raise ValueError("p3点坐标应该是元组、列表或numpy数组类型")

        # 如果点坐标不是np.array类型的, 则将其转化为np.array类型, 否则保持不变
        p0 = np.array(p0) if isinstance(p0, (tuple, list)) else p0
        p1 = np.array(p1) if isinstance(p1, (tuple, list)) else p1
        p2 = np.array(p2) if isinstance(p2, (tuple, list)) else p2
        p3 = np.array(p3) if isinstance(p3, (tuple, list)) else p3

        # 开始插值
        points = []
        for t in np.linspace(0, 1, inserted + 2):
            points.append(p0 * np.power((1 - t), 3) +
                          3 * p1 * t * np.power((1 - t), 2) +
                          3 * p2 * (1 - t) * np.power(t, 2) +
                          p3 * np.power(t, 3))

        return np.array(points)

    @staticmethod
    def smoothing_base_bezier(data_x, data_y, k=0.5, inserted=10, closed=False):
        """基于三阶贝塞尔曲线的数据平滑算法。

        :param data_x: x维度数据集, tuple、list或numpy.ndarray类型
        :param data_y: y维度数据集, tuple、list或numpy.ndarray类型
        :param k: 调整平滑曲线形状的因子, 取值一般在0.2~0.6之间。默认值为0.5
        :param inserted: 两个原始数据点之间插值的数量。默认值为10
        :param closed: 曲线是否封闭, 如是, 则首尾相连。默认曲线不封闭
        :return: 插值完成后的点
        """
        # 如果数据集为list, 则将其转换为np.array
        # 如果本身就是np.array, 则不转换
        # 否则数据集不符合所需数据类型
        if isinstance(data_x, (list, tuple)):
            data_x = np.array(data_x)
        elif not isinstance(data_x, np.ndarray):
            raise ValueError("x数据集应该是元组、列表或numpy数组类型")
        if isinstance(data_y, (list, tuple)):
            data_y = np.array(data_y)
        elif not isinstance(data_y, np.ndarray):
            raise ValueError("y数据集应该是元组、列表或numpy数组类型")

        # 检查数据长度是否一致
        if not data_x.shape == data_y.shape:
            raise ValueError("x数据集和y数据集长度不匹配")

        # 第1步: 生成原始数据折线中点集
        mid_points = []
        for i in range(1, data_x.shape[0]):
            mid_points.append({
                'start': (data_x[i - 1], data_y[i - 1]),
                'end': (data_x[i], data_y[i]),
                'mid': ((data_x[i] + data_x[i - 1]) / 2.0, (data_y[i] + data_y[i - 1]) / 2.0)
            })

        # 当曲线需要闭合时, 如果曲线不需要闭合, 则不会执行
        if closed:
            mid_points.append({
                'start': (data_x[-1], data_y[-1]),
                'end': (data_x[0], data_y[0]),
                'mid': ((data_x[0] + data_x[-1]) / 2.0, (data_y[0] + data_y[-1]) / 2.0)
            })

        # 第2步: 找出中点连线及其分割点
        split_points = []
        for i in range(len(mid_points)):
            if i < (len(mid_points) - 1):
                j = i + 1
            elif closed:
                j = 0
            else:
                continue

            x00, y00 = mid_points[i]['start']
            x01, y01 = mid_points[i]['end']
            x10, y10 = mid_points[j]['start']
            x11, y11 = mid_points[j]['end']
            d0 = np.sqrt(np.power((x00 - x01), 2) + np.power((y00 - y01), 2))
            d1 = np.sqrt(np.power((x10 - x11), 2) + np.power((y10 - y11), 2))
            k_split = 1.0 * d0 / (d0 + d1)

            mx0, my0 = mid_points[i]['mid']
            mx1, my1 = mid_points[j]['mid']

            split_points.append({
                'start': (mx0, my0),
                'end': (mx1, my1),
                'split': (mx0 + (mx1 - mx0) * k_split, my0 + (my1 - my0) * k_split)
            })

        # 第3步: 平移中点连线, 调整端点, 生成控制点
        crt_points = []
        for i in range(len(split_points)):
            vx, vy = mid_points[i]['end']  # 当前顶点的坐标
            dx = vx - split_points[i]['split'][0]  # 平移线段x偏移量
            dy = vy - split_points[i]['split'][1]  # 平移线段y偏移量

            sx, sy = split_points[i]['start'][0] + dx, split_points[i]['start'][1] + dy  # 平移后线段起点坐标
            ex, ey = split_points[i]['end'][0] + dx, split_points[i]['end'][1] + dy  # 平移后线段终点坐标

            cp0 = sx + (vx - sx) * k, sy + (vy - sy) * k  # 控制点坐标
            cp1 = ex + (vx - ex) * k, ey + (vy - ey) * k  # 控制点坐标

            if crt_points:
                crt_points[-1].insert(2, cp0)
            else:
                crt_points.append([mid_points[0]['start'], cp0, mid_points[0]['end']])

            if closed:
                if i < (len(mid_points) - 1):
                    crt_points.append([mid_points[i + 1]['start'], cp1, mid_points[i + 1]['end']])
                else:
                    crt_points[0].insert(1, cp1)
            else:
                if i < (len(mid_points) - 2):
                    crt_points.append([mid_points[i + 1]['start'], cp1, mid_points[i + 1]['end']])
                else:
                    crt_points.append(
                        [mid_points[i + 1]['start'], cp1, mid_points[i + 1]['end'], mid_points[i + 1]['end']])
                    crt_points[0].insert(1, mid_points[0]['start'])

        # 第4步: 应用贝塞尔曲线方程插值
        out = []
        for item in crt_points:
            group = DrawSmoothCurves.bezier_curve(item[0], item[1], item[2], item[3], inserted)
            out.append(group[:-1])

        out.append(group[-1:])
        out = np.vstack(out)

        return out.T[0], out.T[1]

    @staticmethod
    def delete_error_point(filename, x, y):
        """删除错误点。

        :param filename: 文件名, 只有 rocky 文件需要处理, 其他的不需要
        :param x: x轴上所有的点
        :param y: y轴上所有的点
        :return: 处理完成后的x, y点
        """
        need_delete_index_list = []

        if 'rock' in filename:
            for i in range(x.size - 1):
                # 当x值倒退时删除前一个点
                if x[i] < x[i - 1]:
                    need_delete_index_list.append(i - 1)

        return np.delete(x, need_delete_index_list), np.delete(y, need_delete_index_list)


class DataPreProcessing:
    """数据预处理类。"""

    @staticmethod
    def data_preprocessing(df):
        """数据预处理

        :param df: DataFrame类型数据
        :return: 处理完成后的DataFrame
        """
        # u1乘1000, rf1除1000, u1(mm), rf1(KN)
        df['u1(mm)'] = df['u1'] * 1000
        df['rf1(KN)'] = df['rf1'] / 1000

        return df


class YieldDisplacement:
    """计算 屈服位移/屈服点荷载 类"""

    @staticmethod
    def get_intersection(k, x, y, P):
        """获取 kx=P 的 x值在曲线上的坐标。

        :param k: 斜率
        :param x: x轴上所有的点
        :param y: y轴上所有的点
        :param P: P值(P_max 或 P_min)
        :return: 曲线上的坐标
        """
        x_p = P / k  # 计算出 kx=P 时的 x
        x_index = ((x - x_p) >= 0).argmax()  # 计算出比x大的第一个值
        # 该x的值在曲线散点中不存在, 则取x值两侧的点求平均来近似
        if x[x_index] != x_p:
            x_p = np.mean([x[x_index - 1], x[x_index]])
            y_p = np.mean([y[x_index - 1], y[x_index]])
            point = (x_p, y_p)  # 第一个点的坐标
        else:
            point = (x[x_index], y[x_index])

        return point

    @staticmethod
    def get_first_point(x, y, positive=True):
        """获取正方向上的第一个点

        :param x: x轴上所有的点
        :param y: y轴上所有的点
        :param positive: 是否为正方向上的第一个点, 默认为 True
        :return: 正方向上第一个点的坐标
        """
        if positive:
            # 找到第一个可行点(x, y 值都大于 0)
            x_first_index = (x != 0).argmax()  # 找到第一个x值不为0的点
            y_first_index = (y != 0).argmax()  # 找到第一个y值不为0的点
            first_index = max(x_first_index, y_first_index)  # 第一个可行点即为两者坐标较大的那个点
        else:
            # 找到第一个可行点(x, y 值都小于 0)
            x_first_index = (x != 0).argmin() - 1  # 找到第一个x值不为0的点
            y_first_index = (y != 0).argmin() - 1  # 找到第一个y值不为0的点
            first_index = min(x_first_index, y_first_index)  # 第一个可行点即为两者坐标较大的那个点

        # 计算切线斜率
        k = y[first_index] / x[first_index]

        p = y.max() if positive else y.min()
        return YieldDisplacement.get_intersection(k, x, y, p)

    @staticmethod
    def get_smoothing_data(df, filename, inserted=10):
        """获取平滑后的数据

        :param df: DataFrame
        :param filename: 文件名
        :param inserted: 两个原始数据点之间插值的数量。默认值为10
        :return: 平滑后的数据点
        """
        # 提取数据集 x, y
        x, y = df['u1(mm)'].values, df['rf1(KN)'].values
        x, y = DrawSmoothCurves.delete_error_point(filename, x, y)

        # 利用 三阶贝塞尔曲线的数据平滑算法 对曲线进行平滑处理
        x_curve, y_curve = DrawSmoothCurves.smoothing_base_bezier(x, y, k=0.5, inserted=inserted, closed=False)

        return x_curve, y_curve

    @staticmethod
    def get_yield_displacement(df, filename):
        """计算 屈服位移/屈服点荷载

        :param df: DataFrame类型数据
        :param filename: 文件名
        :return: 处理完成后的DataFrame
        """
        x, y = YieldDisplacement.get_smoothing_data(df, filename)
        df['屈服位移'] = None
        df['屈服点荷载'] = None

        # x轴正方向的点坐标
        positive_index = np.where(x >= 0)
        x_curve_positive = x[positive_index]
        y_curve_positive = y[positive_index]

        # 获得第一个点
        first_positive_point = YieldDisplacement.get_first_point(x_curve_positive, y_curve_positive, True)
        # 根据第一个点获得第二个点
        k = first_positive_point[1] / first_positive_point[0]
        P_y = YieldDisplacement.get_intersection(k, x_curve_positive, y_curve_positive, y_curve_positive.max())
        df['屈服位移'].iloc[0] = P_y[0]
        df['屈服点荷载'].iloc[0] = P_y[1]

        # **********************************************************************
        # x轴负方向的点坐标
        positive_index = np.where(x <= 0)
        x_curve_positive = x[positive_index]
        y_curve_positive = y[positive_index]

        # 获得第一个点
        first_positive_point = YieldDisplacement.get_first_point(x_curve_positive, y_curve_positive, False)
        # 根据第一个点获得第二个点
        k = first_positive_point[1] / first_positive_point[0]
        P_y = YieldDisplacement.get_intersection(k, x_curve_positive, y_curve_positive, y_curve_positive.min())
        df['屈服位移'].iloc[1] = P_y[0]
        df['屈服点荷载'].iloc[1] = P_y[1]

        return df


class GetMaxMin:
    """计算 极限位移/极限荷载"""

    @staticmethod
    def get_max_min(df, df_rock):
        """计算 极限位移/极限荷载

        :param df: DataFrame类型数据
        :param df_rock: DataFrame类型数据
        :return: 处理完成后的df_rock
        """
        # 创建 极限位移 与 极限荷载 列, 默认赋值为 空
        df_rock['极限位移'] = None
        df_rock['极限荷载'] = None

        # 将每一列的第一个元素填充为所需元素
        # df['u1'].iloc[i] 指取 u1 列第 i + 1 个元素(因为 i 从 0 开始)
        # argmax() 可以获得最大值元素的下标
        df_rock['极限位移'].iloc[0] = df['u1(mm)'].iloc[df['u1(mm)'].abs().argmax()]
        df_rock['极限荷载'].iloc[0] = df['rf1(KN)'].iloc[df['rf1(KN)'].abs().argmax()]

        return df_rock


class DisplacementDuctilityFactor:
    """计算 位移延性系数"""

    @staticmethod
    def calc_displacement_ductility_factor(df_rock):
        """位移延性系数。

        :param df_rock: DataFrame类型数据
        :return: 处理完成后的df_rock
        """
        df_rock['位移延性系数'] = (df_rock['极限位移'].iloc[0] / df_rock['屈服位移']).abs()

        return df_rock


class EnergyDissipationCoefficient:
    """计算 能量耗散系数"""

    @staticmethod
    def calc_energy_dissipation_coefficient(x_curves, y_curves, show_figure=False):
        """计算 能量耗散系数。

        :param x_curves: x轴上所有的点
        :param y_curves: y轴上所有的点
        :param show_figure: 是否显示曲形图像, 默认不显示
        :return: 能量耗散系数 E
        """
        # 头尾连结, 使之成为一个封闭区域
        x_curves = np.r_[x_curves, x_curves[0]]
        y_curves = np.r_[y_curves, y_curves[0]]

        if show_figure:
            plt.plot(x_curves, y_curves)
            plt.show()

        # 找出两个三角形的点
        x_max_index, x_min_index = x_curves.argmax(), x_curves.argmin()
        x_max, x_min = x_curves[x_max_index], x_curves[x_min_index]
        y_max, y_min = y_curves[x_max_index], y_curves[x_min_index]

        # 计算面积
        S_area = np.abs(np.trapz(y_curves, x_curves))
        S_triangle = np.abs((x_max * y_max + x_min * y_min) / 2)
        E = S_area / S_triangle

        return E

    @staticmethod
    def get_index_list(x):
        """找到可以作为终点的点。"""
        index_list = []
        for i in range((x > 0).argmax(), x.size - 1):
            if x[i] < 0 <= x[i + 1]:
                # 第i + 1个为止可以作为起点或终点
                index_list.append(i + 1)

        return index_list


class Main:
    def __init__(self):
        """初始化类"""
        # .csv文件输入目录
        self.input_dir = './result'
        # 读取 data 目录下所有 .csv 文件
        self.total_filename_list = os.listdir(self.input_dir)
        self.rock_filename_list = [filename for filename in self.total_filename_list if filename.endswith("rock.csv")]

        # csv输出目录
        self.csv_output_dir = './result'
        # 如果当前目录下不存在 output 目录, 则创建该目录
        if not os.path.exists(self.csv_output_dir):
            os.mkdir(self.csv_output_dir)

        # 图片输出目录
        self.png_output_dir = './output'
        # 如果当前目录下不存在 output 目录, 则创建该目录
        if not os.path.exists(self.png_output_dir):
            os.mkdir(self.png_output_dir)

    def get_df(self, filename):
        """根据 filename 读取 csv 文件

        :param filename: 文件名
        :return: DataFrame
        """
        # 获取 .csv 文件地址
        file_path = os.path.join(self.input_dir, filename)

        return pd.read_csv(file_path, encoding='gbk')

    def save_df(self, df, filename):
        """根据 filename 保存 csv 文件

        :param df: DataFrame
        :param filename: 文件名
        """
        # 获取输出文件地址
        output_file_path = os.path.join(self.csv_output_dir, filename)
        df.to_csv(output_file_path, index=False, encoding='gbk')

    def data_preprocessing(self):
        """数据预处理"""
        # 遍历所有 .csv 文件
        for filename in tqdm(self.total_filename_list, desc="数据预处理中", leave=False):
            df = self.get_df(filename)
            df = DataPreProcessing.data_preprocessing(df)
            self.save_df(df, filename)

    def draw_smooth_curves(self):
        """绘制平滑曲线"""
        matplotlib.rcParams['font.family'] = 'sans-serif'  # 修改字体
        matplotlib.rcParams['font.sans-serif'] = 'Times New Roman, NSimSun'  # 修改字体
        figure_size = (8, 5)    # 宽度，高度，单位英寸
        title_font_size = 16    # 标题字体大小
        line_width = 1          # 线条宽度
        color = 'blue'          # 线条颜色, 所有颜色请详见 https://matplotlib.org/examples/color/named_colors.html ,或使用形如 #000000 方式调色
        line_style = 'solid'    # 线条类型, 所有类型请详见 https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/linestyles.html
        x_label = "∆/mm"        # x轴标注
        y_label = "F/KN"        # y轴标注
        x_label_font_size = 16  # x轴标注字体大小
        y_label_font_size = 16  # y轴标注字体大小
        x_ticks_font_size = 16  # x轴刻度字体大小
        y_ticks_font_size = 16  # y轴刻度字体大小
        legend = "rf1(KN)"      # 图例
        legend_font_size = 16   # 图例字体大小

        # 遍历所有 .csv 文件
        for filename in tqdm(self.total_filename_list, desc="绘制平滑曲线中", leave=False):
            df = self.get_df(filename)

            # 提取数据集 x, y
            x, y = df['u1(mm)'].values, df['rf1(KN)'].values
            x, y = DrawSmoothCurves.delete_error_point(filename, x, y)
            # 利用 三阶贝塞尔曲线的数据平滑算法 对曲线进行平滑处理
            x_curve, y_curve = DrawSmoothCurves.smoothing_base_bezier(x, y, k=0.5, inserted=10, closed=False)

            x_limits = (x.min() * 1.1, x.max() * 1.1)  # x轴范围
            y_limits = (y.min() * 1.1, y.max() * 1.1)  # y轴范围
            # 开始一个新的图像
            plt.figure(figsize=figure_size)
            # plt.title("rf1(KN)", fontsize=title_font_size)
            # 绘制图像
            plt.rcParams['xtick.direction'] = 'in'#将x轴的刻度线方向设置向内
            plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度线方向设置向内
            plt.plot(x_curve, y_curve, linewidth=line_width, linestyle=line_style, color=color)
            plt.xticks(fontsize=x_ticks_font_size)
            plt.yticks(fontsize=y_ticks_font_size)
            plt.xlim(x_limits)
            plt.ylim(y_limits)
            if filename.endswith('rock.csv'):
                ax = plt.gca()
                ax.spines['top'].set_color('none')
                ax.spines['right'].set_color('none')
                ax.xaxis.set_ticks_position('bottom')
                ax.spines['bottom'].set_position(('data', 0))
                ax.yaxis.set_ticks_position('left')
                ax.spines['left'].set_position(('data', 0))
                plt.xticks([-5, -3, -1, 0, 1, 3, 5])
                plt.yticks([-120, -90, -60, -30, 30, 60, 90, 120])
                plt.xlabel(x_label, fontsize=x_label_font_size, horizontalalignment='right', position=(1, 1))
                plt.ylabel(y_label, fontsize=y_label_font_size, verticalalignment='top', position=(1, 1))
                plt.legend((legend,), fontsize=legend_font_size)
            else:
                plt.xticks([-5, -3, -1, 1, 3, 5])
                plt.yticks([-120, -90, -60, -30, 0, 30, 60, 90, 120])
                plt.xlabel(x_label, fontsize=x_label_font_size)
                plt.ylabel(y_label, fontsize=y_label_font_size)
                plt.grid(axis='x', linestyle='-')
                plt.grid(axis='y', linestyle='-')
            # 获取输出文件地址
            output_file_path = os.path.join(self.png_output_dir, filename.replace("csv", "png"))
            # 保存图像
            plt.savefig(output_file_path)

    def calc_yield_displacement(self):
        """计算 屈服位移"""
        # 遍历所有 .csv 文件
        for filename in tqdm(self.rock_filename_list, desc="计算 屈服位移 中", leave=False):
            df = self.get_df(filename)

            df = YieldDisplacement.get_yield_displacement(df, filename)

            self.save_df(df, filename)

    def calc_max_min_and_displacement_ductility_factor(self):
        """计算 极限位移/极限荷载/位移延性系数"""
        # 遍历所有 .csv 文件
        for rock_filename in tqdm(self.rock_filename_list, desc="计算 极限位移/极限荷载/位移延性系数 中", leave=False):
            filename = rock_filename.replace('-rock', '')
            df = self.get_df(filename)
            df_rock = self.get_df(rock_filename)

            df_rock = GetMaxMin.get_max_min(df, df_rock)
            df_rock = DisplacementDuctilityFactor.calc_displacement_ductility_factor(df_rock)

            self.save_df(df_rock, rock_filename)

    def calc_energy_dissipation_coefficient(self):
        """计算 能量耗散系数。"""
        # 遍历所有 .csv 文件
        for filename in tqdm(self.rock_filename_list, desc="计算 能量耗散系数中", leave=False):
            without_rock_filename = filename.replace('-rock', '')
            df = self.get_df(without_rock_filename)
            df_rock = self.get_df(filename)
            df_rock['能量耗散系数'] = None
            df_rock['等效粘滞阻尼系数'] = None

            # 平滑后数据
            x, y = YieldDisplacement.get_smoothing_data(df, without_rock_filename, inserted=10)
            # 原始数据, 直接算
            # x, y = df['u1(mm)'].to_numpy(), df['rf1(KN)'].to_numpy()

            index_list = EnergyDissipationCoefficient.get_index_list(x)

            E_list = []
            start_index = 0
            for end_index in index_list:
                # 获取曲形的点
                x_curves = x[start_index: end_index + 1]
                y_curves = y[start_index: end_index + 1]
                # 下个曲形的起点为当前曲形的终点
                start_index = end_index

                # 计算 E
                E = EnergyDissipationCoefficient.calc_energy_dissipation_coefficient(x_curves, y_curves)
                E_list.append(E)

            # start_index指向的不是最后一个数据, 说明x中还存在数据
            if start_index != (x.size - 1):
                x_curves = x[start_index:]
                y_curves = y[start_index:]
                last_index_list = np.where(x_curves >= 0)

                x_curves = x_curves[last_index_list]
                y_curves = y_curves[last_index_list]

                E = EnergyDissipationCoefficient.calc_energy_dissipation_coefficient(x_curves, y_curves)
                E_list.append(E)

            df_rock['能量耗散系数'].iloc[0: len(E_list)] = E_list
            df_rock['等效粘滞阻尼系数'] = df_rock['能量耗散系数'] / (2 * np.pi)

            self.save_df(df_rock, filename)


class Pipeline:
    @staticmethod
    def begin_tasks():
        sleep_time = 0.5
        main = Main()
        task_count = len(dir(main)) - 33

        pbar = tqdm(total=task_count, desc="开始处理任务")
        # 数据预处理
        main.data_preprocessing()
        pbar.update(1)
        time.sleep(sleep_time)

        # 绘制平滑曲线
        main.draw_smooth_curves()
        pbar.update(1)
        time.sleep(sleep_time)

        # 绘制平滑曲线
        main.calc_yield_displacement()
        pbar.update(1)
        time.sleep(sleep_time)

        # 计算 极限位移/极限荷载/位移延性系数
        main.calc_max_min_and_displacement_ductility_factor()
        pbar.update(1)
        time.sleep(sleep_time)

        # 计算 能量耗散系数
        main.calc_energy_dissipation_coefficient()
        pbar.update(1)
        time.sleep(sleep_time)

        pbar.close()
        time.sleep(sleep_time)
        print("所有任务全部完成")


if __name__ == '__main__':
    import warnings

    # 忽略警告
    warnings.filterwarnings("ignore")
    # 开始执行任务
    Pipeline.begin_tasks()
