#!/usr/bin/env bash

set -u

if [ "$#" -ne 1 ]; then
	echo "usage: $0 <tasks.txt>" >&2
	exit 2
fi

TASKS=$1
BIN=$(dirname "$0")/md5fastcoll

if [ ! -r "$TASKS" ]; then
	echo "$0: cannot read task file: $TASKS" >&2
	exit 1
fi

if [ ! -x "$BIN" ]; then
	echo "$0: $BIN not found or not executable; run 'make all' first" >&2
	exit 1
fi

status=0
n=0

# Fields are space-separated: <input_file> <output_file1> <output_file2>.
# 主路径：单进程多线程批量模式
if "$BIN" --tasks "$TASKS" >/dev/null; then
	n=$(awk 'NF { c++ } END { print c+0 }' "$TASKS")    # 非空行数
	status=0
else
	# 降级路径：主路径整体失败就退回原来的逐个串行方式
	echo "$0: bulk mode failed, falling back to serial" >&2
	status=0
	while read -r infile out1 out2 rest; do
		[ -n "${infile:-}" ] || continue
		if ! "$BIN" -q -p "$infile" -o "$out1" "$out2" >/dev/null; then
			echo "$0: failed on $infile" >&2
			status=1
		else
			n=$((n + 1))
		fi
	done < "$TASKS"
fi


echo "$0: generated $n collisions"
exit "$status"
