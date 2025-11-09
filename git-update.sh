#!/bin/bash

# ==============================================================================
# Git 快速更新脚本 (git-update.sh)
# 用法: 在你的 Git 仓库根目录下运行 ./git-update.sh
# 功能: 自动执行 add, commit, 和 push 的标准流程。
# ==============================================================================

# --- 颜色定义 (让输出更好看) ---
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# --- 步骤 1: 检查是否在 Git 仓库中 ---
if ! git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
    echo -e "${RED}错误: 当前目录不是一个 Git 仓库。${NC}"
    exit 1
fi

echo -e "${GREEN}=== Git 快速更新脚本 ===${NC}"

# --- 步骤 2: 显示当前状态 ---
echo -e "\n${YELLOW}--- 1. 当前仓库状态 ---${NC}"
git status
echo -e "${YELLOW}--------------------------${NC}"

# --- 确认是否继续 ---
read -p "以上是当前的更改，是否继续？(y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    echo -e "${RED}操作已取消。${NC}"
    exit 1
fi

# --- 步骤 3: 添加所有文件 ---
echo -e "\n${YELLOW}--- 2. 添加所有更改到暂存区 (git add .) ---${NC}"
git add .
echo -e "${GREEN}✅ 所有更改已添加。${NC}"

# --- 步骤 4 & 5: 打开编辑器进行提交 ---
echo -e "\n${YELLOW}--- 3. 即将打开文本编辑器以输入提交信息 ---${NC}"
echo "提示: "
echo "  - 在编辑器中详细输入您的提交信息。"
echo "  - 以 '#' 开头的行是注释，将被忽略。"
echo "  - 保存并关闭编辑器以完成提交。不输入任何信息直接关闭则会取消提交。"
read -p "按回车键继续..."

git commit

# 检查 commit 是否成功
if [ $? -ne 0 ]; then
    echo -e "${RED}错误: Git commit 失败。请检查错误信息。${NC}"
    exit 1
fi

# --- 步骤 6: 推送 ---
echo -e "\n${YELLOW}--- 5. 推送到远程仓库 (git push) ---${NC}"
git push

# 检查 push 是否成功
if [ $? -eq 0 ]; then
    echo -e "\n${GREEN}🚀 更新成功！所有更改已推送到远程仓库。${NC}"
else
    echo -e "\n${RED}错误: Git push 失败。请检查网络连接或权限问题。${NC}"
fi

exit 0
