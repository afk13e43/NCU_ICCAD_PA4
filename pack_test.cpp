#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

// --- 1. 資料結構定義 ---

// 模組結構：儲存尺寸和最終計算出的座標
struct Module
{
    std::string name;
    double w, h; // 寬、高
    double x, y; // 左下角座標 (由 Packer 計算)

    Module(std::string n, double width, double height)
        : name(n), w(width), h(height), x(0), y(0) {}
};

// B*-Tree 節點結構
struct BTreeNode
{
    Module *module;
    BTreeNode *left, *right, *parent;

    BTreeNode(Module *m)
        : module(m), left(nullptr), right(nullptr), parent(nullptr) {}
};

// 輪廓線節點結構
struct ContourNode
{
    double p1, p2; // p1: 輪廓線段的起始點, p2: 該線段的高度/位置
    ContourNode *prev, *next;
};

// --- 2. 輪廓線類別 ---

// 水平輪廓線 (天際線)，記錄 y = f(x)
class HorizontalContour
{
public:
    HorizontalContour()
    {
        // 初始輪廓線是 y=0 的地板
        head = new ContourNode{0.0, 0.0, nullptr, nullptr};
    }

    ~HorizontalContour()
    {
        // ... 釋放鏈結串列記憶體 ...
    }

    // 查詢 [x_start, x_end) 區間內的最大 y 高度
    double get_max_y(double x_start, double x_end)
    {
        if (x_start >= x_end)
            return 0.0;
        double max_y = 0.0;
        ContourNode *it = head;
        // 找到包含 x_start 的線段
        while (it->next != nullptr && it->next->p1 <= x_start)
        {
            it = it->next;
        }
        // 遍歷所有重疊的線段
        while (it != nullptr && it->p1 < x_end)
        {
            max_y = std::max(max_y, it->p2);
            it = it->next;
        }
        return max_y;
    }

    // 用新模組的頂部邊界來更新輪廓線
    void update(double x_start, double x_end, double new_y)
    {
        double original_y_at_x_end = get_max_y(x_end, x_end + 1e-9);

        // 找到插入/修改的起始位置
        ContourNode *it = head;
        while (it->next != nullptr && it->next->p1 < x_start)
        {
            it = it->next;
        }

        // 刪除被完全覆蓋的舊節點
        ContourNode *current = it->next;
        while (current != nullptr && current->p1 < x_end)
        {
            ContourNode *to_delete = current;
            current = current->next;
            // 從鏈結串列中移除
            to_delete->prev->next = to_delete->next;
            if (to_delete->next)
                to_delete->next->prev = to_delete->prev;
            delete to_delete;
        }

        // 插入新的起始節點
        if (it->p1 == x_start)
        { // 如果x座標重合，直接更新
            it->p2 = std::max(it->p2, new_y);
        }
        else
        { // 否則插入新節點
            ContourNode *new_node = new ContourNode{x_start, new_y, it, it->next};
            if (it->next)
                it->next->prev = new_node;
            it->next = new_node;
            it = new_node;
        }

        // 插入新的結束節點以恢復原始高度
        if (!it->next || it->next->p1 > x_end)
        {
            ContourNode *end_node = new ContourNode{x_end, original_y_at_x_end, it, it->next};
            if (it->next)
                it->next->prev = end_node;
            it->next = end_node;
        }
    }

private:
    ContourNode *head;
};

// 垂直輪廓線 (右側邊界)，記錄 x = f(y)
// 邏輯與 HorizontalContour 完全對稱
class VerticalContour
{
public:
    VerticalContour()
    {
        // 初始輪廓線是 x=0 的牆面
        head = new ContourNode{0.0, 0.0, nullptr, nullptr};
    }
    ~VerticalContour()
    {
        // ... 釋放記憶體 ...
    }

    // 查詢 [y_start, y_end) 區間內的最大 x 位置
    double get_max_x(double y_start, double y_end)
    {
        if (y_start >= y_end)
            return 0.0;
        double max_x = 0.0;
        ContourNode *it = head;
        while (it->next != nullptr && it->next->p1 <= y_start)
        {
            it = it->next;
        }
        while (it != nullptr && it->p1 < y_end)
        {
            max_x = std::max(max_x, it->p2);
            it = it->next;
        }
        return max_x;
    }

    // 用新模組的右側邊界來更新輪廓線
    void update(double y_start, double y_end, double new_x)
    {
        // 邏輯與 HorizontalContour::update 對稱
        double original_x_at_y_end = get_max_x(y_end, y_end + 1e-9);
        ContourNode *it = head;
        while (it->next != nullptr && it->next->p1 < y_start)
        {
            it = it->next;
        }
        ContourNode *current = it->next;
        while (current != nullptr && current->p1 < y_end)
        {
            ContourNode *to_delete = current;
            current = current->next;
            to_delete->prev->next = to_delete->next;
            if (to_delete->next)
                to_delete->next->prev = to_delete->prev;
            delete to_delete;
        }
        if (it->p1 == y_start)
        {
            it->p2 = std::max(it->p2, new_x);
        }
        else
        {
            ContourNode *new_node = new ContourNode{y_start, new_x, it, it->next};
            if (it->next)
                it->next->prev = new_node;
            it->next = new_node;
            it = new_node;
        }
        if (!it->next || it->next->p1 > y_end)
        {
            ContourNode *end_node = new ContourNode{y_end, original_x_at_y_end, it, it->next};
            if (it->next)
                it->next->prev = end_node;
            it->next = end_node;
        }
    }

private:
    ContourNode *head;
};

// --- 3. 打包器類別 ---

class BStarTreePacker
{
public:
    void pack(BTreeNode *root)
    {
        if (!root)
            return;

        h_contour = new HorizontalContour();
        v_contour = new VerticalContour();

        // 從根節點開始遞迴打包
        recursive_pack(root);

        delete h_contour;
        delete v_contour;
    }

private:
    HorizontalContour *h_contour;
    VerticalContour *v_contour;

    void recursive_pack(BTreeNode *node)
    {
        if (node == nullptr)
            return;

        Module *m = node->module;
        BTreeNode *p = node->parent;

        double final_x = 0;
        double final_y = 0;

        if (p == nullptr)
        { // --- Case 1: 根節點 ---
            final_x = 0;
            final_y = 0;
        }
        else if (node == p->left)
        { // --- Case 2: 左子節點 (置於父節點之上) ---
            // x 座標與父節點對齊
            final_x = p->module->x;
            // y 座標是父節點上方，以及該 x 區間內其他模組的最高點
            final_y = h_contour->get_max_y(final_x, final_x + m->w);
        }
        else
        { // --- Case 3: 右子節點 (置於父節點之右) ---
            // x 座標是父節點右方，以及該 y 區間內其他模組的最右點
            final_x = v_contour->get_max_x(0, m->h);
            final_x = std::max(final_x, p->module->x + p->module->w); // 必須在父節點右邊

            // y 座標是該 x 區間的最高點 (即盡量靠地)
            final_y = h_contour->get_max_y(final_x, final_x + m->w);
        }

        // 設定最終座標
        m->x = final_x;
        m->y = final_y;

        // 更新兩個全局輪廓線
        h_contour->update(m->x, m->x + m->w, m->y + m->h);
        v_contour->update(m->y, m->y + m->h, m->x + m->w);

        // 遞迴處理子節點
        recursive_pack(node->left);
        recursive_pack(node->right);
    }
};

// --- 4. 主函式：建立與測試 ---

int main()
{
    // 建立模組
    std::map<std::string, Module> modules;
    modules.emplace("A", Module("A", 10, 20));
    modules.emplace("B", Module("B", 8, 15));
    modules.emplace("C", Module("C", 12, 18));
    modules.emplace("D", Module("D", 7, 10));

    // 建立 B*-Tree
    //      A
    //     / \
    //    B   C
    //       /
    //      D
    BTreeNode *root = new BTreeNode(&modules.at("A"));
    root->left = new BTreeNode(&modules.at("B"));
    root->left->parent = root;

    root->right = new BTreeNode(&modules.at("C"));
    root->right->parent = root;

    root->right->left = new BTreeNode(&modules.at("D"));
    root->right->left->parent = root->right;

    // 建立並執行打包器
    BStarTreePacker packer;
    packer.pack(root);

    // 輸出結果
    std::cout << "Final Placements:" << std::endl;
    for (auto const &[name, module] : modules)
    {
        std::cout << "Module " << name
                  << " (w:" << module.w << ", h:" << module.h << ")"
                  << " -> placed at (" << module.x << ", " << module.y << ")" << std::endl;
    }

    // 清理 B*-Tree 記憶體 (在真實應用中需要完整的清理函式)
    delete root->right->left;
    delete root->right;
    delete root->left;
    delete root;

    return 0;
}