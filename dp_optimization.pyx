import numpy as np
cimport numpy as np
from cpython.pycapsule cimport PyObject  # 导入 PyObject 类型用于字符串数组

# 修改类型声明
def run_dp(np.ndarray[np.float64_t] filter_seq, np.ndarray[np.int64_t] merged_df_start, np.ndarray[np.int64_t] merged_df_end,
           np.ndarray[np.int64_t] merged_df_distance, np.ndarray[np.object_] merged_df_motif, int length):

    # 定义常量
    cdef double gap_penalty = 1
    cdef double distance_penalty = 1
    cdef double perfect_bonus = 0.5

    # 初始化 dp 数组和 pre 数组
    cdef np.ndarray[np.float64_t] dp = np.zeros(length + 1, dtype=np.float64)
    cdef np.ndarray[np.int64_t] pre = np.full((length + 1, 3), None, dtype=np.int64)  # pre_i, motif_id, motif
    cdef int i, idx = 0
    cdef int pre_i
    cdef double motif_len, distance, bonus

    # 将 merged_df_motif 转换为列表
    cdef list motif_list = [merged_df_motif[i] for i in range(len(merged_df_motif))]

    for i in range(1, length + 1):
        # 跳过一个基因
        if dp[i-1] - gap_penalty > 0:
            dp[i] = dp[i-1] - gap_penalty
            pre[i] = (i-1, None, None)

        while idx < len(merged_df_start) and merged_df_end[idx] <= i:
            if merged_df_end[idx] < i:
                idx += 1
                continue

            pre_i = merged_df_start[idx]
            motif_len = len(motif_list[idx])  # 计算 motif 的长度
            distance = merged_df_distance[idx]
            bonus = perfect_bonus * motif_len if distance == 0 else 0

            if dp[pre_i] + motif_len - distance * distance_penalty + bonus >= dp[i]:
                dp[i] = dp[pre_i] + motif_len - distance * distance_penalty + bonus
                pre[i] = (pre_i, idx, motif_list[idx])

            idx += 1

        if i % 5000 == 0:
            print(f'DP: {i // 1000} kbp is Done!')

    print('DP complete!')
