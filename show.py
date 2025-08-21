import matplotlib.pyplot as plt
import matplotlib.patches as patches

# === 1. 讀檔 ===
RESULT_FILE = "result1000.txt"
detail=False
with open(RESULT_FILE, "r") as f:
    lines = [l.strip() for l in f if l.strip()]

chip_area   = float(lines[0])              # 第 1 行
chip_w, chip_h = map(float, lines[1].split())  # 第 2 行
inl_value   = float(lines[2])              # 第 3 行

block_lines = lines[3:]                    # 其餘各行

blocks = []
for l in block_lines:
    # 範例：MM0 0.00 0.00 (9.11 5.54 8 1)
    name, x, y, rest = l.split(maxsplit=3)
    x, y = float(x), float(y)
    w, h, cm, rm = map(float, rest.strip("()").split()[:4])
    blocks.append((name, x, y, w, h))

# === 2. 畫圖 ===
#plt.figure(figsize=(chip_w/4, chip_h/4))   # 對比例縮小，避免圖太大
fig, ax = plt.subplots(figsize=(10, 18)) 
ax = plt.gca()
ax.set_aspect('equal')
ax.set_xlim(0, chip_w)
ax.set_ylim(0, chip_h)

# 晶片外框
ax.add_patch(
    patches.Rectangle((0, 0), chip_w, chip_h,
                      fill=False, edgecolor='black',
                      linewidth=2, linestyle='--',
                      label=f"Chip {chip_w:.2f}×{chip_h:.2f}")
)

# 每個 block
for name, x, y, w, h in blocks:
    rect = patches.Rectangle((x, y), w, h,
                             linewidth=1, edgecolor='tab:blue',
                             facecolor='skyblue', alpha=0.5)
    ax.add_patch(rect)
    if detail:
        ax.text(x + w/2, y + h/2,
            f"{name}\n{w:.1f}×{h:.1f}",
            ha='center', va='center', fontsize=8)
    

# 標題／格線
plt.title(f"PA4 Floorplan  /  Cost={chip_area:.2f}  /  INL={inl_value:.2f}")
plt.grid(True, linestyle=":", linewidth=0.5)

plt.show()
