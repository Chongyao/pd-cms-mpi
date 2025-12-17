#!/bin/bash

# 尝试自动获取当前进程的 MPI Rank ID
# OpenMPI 使用 OMPI_COMM_WORLD_RANK
# Intel MPI / MPICH 使用 PMI_RANK
RANK="0"
if [ -n "$OMPI_COMM_WORLD_RANK" ]; then
    RANK=$OMPI_COMM_WORLD_RANK
elif [ -n "$PMI_RANK" ]; then
    RANK=$PMI_RANK
elif [ -n "$SLURM_PROCID" ]; then
    RANK=$SLURM_PROCID
fi

# 定义输出文件名，例如: time_log_0.txt, time_log_1.txt ...
LOG_FILE="time_log_${RANK}.txt"

# 执行命令，并将 time 的输出 (-o) 重定向到文件
# "$@" 代表传递给脚本的所有参数（即你的程序和参数）
/usr/bin/time -v -o "$LOG_FILE" "$@"