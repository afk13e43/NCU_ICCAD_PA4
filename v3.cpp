#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <queue>
#include <set>
#include <cmath>
#include <limits>
#include <random> // HINT: 為了模擬退火中的隨機過程，需要這個標頭檔
#include <list>

using namespace std;
// 用於儲存單一形態的幾何資訊
// 格式: (width height col_multiple row_multiple)
struct Shape
{
    double width;
    double height;
    int col_multiple;
    int row_multiple;
};

// 用於儲存一個裝置 (模組) 的資訊
// 包含裝置名稱和一個或多個可能的形態
struct Module
{
    std::string name;
    std::vector<Shape> possible_shapes;
};

class Placer
{
public:
    // --- 需要實作 (1): 建構子 ---
    Placer(const std::vector<Module> &modules)
    {
        // 初始化 root 指標
        root = nullptr;
        Node *prev_node = nullptr;

        // 1. 遍歷傳入的 modules vector
        // 1. 建立所有節點物件，但不建立樹連接
        for (const auto &module : modules)
        {
            // 在 Placer 建構子中
            const Shape &initial_shape = module.possible_shapes[0];
            // [修正] 乘上倍率來計算總尺寸
            double total_width = initial_shape.width * initial_shape.col_multiple;
            double total_height = initial_shape.height * initial_shape.row_multiple;
            Node *new_node = new Node(module.name, total_width, total_height);
            new_node->module_info = &module;
            this->nodes.push_back(new_node);
        }
        // 2. 隨機建立 B*-Tree 的拓撲結構
        root = nodes[0]; // 第一個節點總是根

        // 從第二個節點開始，為每個節點隨機找一個父節點
        for (size_t i = 1; i < nodes.size(); ++i)
        {
            Node *new_node = nodes[i];
            bool placed = false;

            // 持續嘗試，直到為 new_node 找到一個家
            while (!placed)
            {
                // 從已經在樹中的節點 (0 到 i-1) 中隨機選一個當爸爸
                int parent_idx = rand() % i;
                Node *parent_node = nodes[parent_idx];

                // 隨機決定要放左邊還是右邊
                if (rand() % 2 == 0)
                {
                    // 嘗試放在左邊
                    if (parent_node->left_child == nullptr)
                    {
                        parent_node->left_child = new_node;
                        new_node->parent = parent_node;
                        placed = true;
                    }
                    // 如果左邊被佔了，試試右邊
                    else if (parent_node->right_child == nullptr)
                    {
                        parent_node->right_child = new_node;
                        new_node->parent = parent_node;
                        placed = true;
                    }
                }
                else
                {
                    // 嘗試放在右邊
                    if (parent_node->right_child == nullptr)
                    {
                        parent_node->right_child = new_node;
                        new_node->parent = parent_node;
                        placed = true;
                    }
                    // 如果右邊被佔了，試試左邊
                    else if (parent_node->left_child == nullptr)
                    {
                        parent_node->left_child = new_node;
                        new_node->parent = parent_node;
                        placed = true;
                    }
                }
                // 如果選中的父節點左右都滿了，迴圈會繼續，重新選另一個父節點
            }
        }
        // 初始化最佳解的儲存空間 (可以先不做，等 SA 寫好再加)
        best_root = nullptr;
    }

    void solve()
    {
        // 1. 初始參數設定
        double T = 10000.0;               // 初始溫度
        double T_min = 0.1;               // 終止溫度
        double cooling_rate = 0.99;       // 降溫速率
        int max_iter = nodes.size() * 10; // 每個溫度下的迭代次數

        // 2. 產生初始解並計算成本
        packing();
        double current_cost = calculate_cost();
        double best_cost = current_cost;
        copy_tree(nodes, root, best_nodes_state, best_root); // 備份初始最佳解

        // 3. 模擬退火主迴圈
        while (T > T_min)
        {
            for (int i = 0; i < max_iter; ++i)
            {
                // a. 儲存當前狀態 (以便退回)
                std::vector<Node *> backup_nodes_state;
                Node *backup_root = root;
                copy_tree(nodes, root, backup_nodes_state, backup_root);

                // b. 擾動 (Perturb)，產生新解
                int move_type = rand() % 3; // 隨機選擇一種擾動
                perturb(move_type);

                // c. 計算新解的成本
                packing(); // 擾動後需要重新 packing
                double new_cost = calculate_cost();

                // d. 判斷是否接受新解
                double delta_cost = new_cost - current_cost;
                if (delta_cost < 0)
                { // 新解更好，直接接受
                    current_cost = new_cost;
                    if (new_cost < best_cost)
                    { // 更新全域最佳解
                        best_cost = new_cost;
                        copy_tree(nodes, root, best_nodes_state, best_root);
                    }
                }
                else
                { // 新解比較差，以一定機率接受
                    double accept_prob = exp(-delta_cost / T);
                    if ((double)rand() / RAND_MAX < accept_prob)
                    {
                        current_cost = new_cost;
                    }
                    else
                    {
                        // 不接受，還原到擾動前的狀態
                        restore_tree(backup_nodes_state, backup_root, nodes, root);
                    }
                }
            }
            // e. 降溫
            cout << "Current temperature: " << T << ", Best cost: " << best_cost << endl;
            T *= cooling_rate;
        }

        // 4. 迴圈結束，還原到找到的最佳解
        restore_tree(best_nodes_state, best_root, nodes, root);
        packing(); // 確保最終輸出的是最佳解的 packing 結果
        cout << "SA finished. Best cost found: " << best_cost << endl;
    }

    // --- 已為您實作 (3): 檔案輸出 ---
    // --- 已為您實作 (3): 檔案輸出 (修正版) ---
    void write_output(const std::string &filename)
    {
        // 1. 計算晶片尺寸
        double chip_width = this->current_chip_width;
        double chip_height = this->current_chip_height;
        double chip_area = chip_width * chip_height;

        // 2. 計算 INL
        double inl = calculate_inl();

        // 3. 開啟 ofstream
        ofstream out_file(filename);
        if (!out_file.is_open())
        {
            cerr << "錯誤: 無法開啟輸出檔案 " << filename << endl;
            return;
        }

        // 4. 寫入前三行 (面積、尺寸、INL)，這部分不變
        out_file << fixed << setprecision(4) << chip_area << endl;
        out_file << fixed << setprecision(2) << chip_width << " " << chip_height << endl;
        out_file << fixed << setprecision(2) << inl << endl;

        // --- 【新增修正】 ---
        // 5. 為了輸出排序好的模組列表，建立一個 nodes 的副本
        std::vector<Node *> sorted_output_nodes = this->nodes;

        // 6. 使用與 INL 計算中完全相同的 "自然排序" 邏輯對副本進行排序
        sort(sorted_output_nodes.begin(), sorted_output_nodes.end(), [](const Node *a, const Node *b)
             {
        size_t first_digit_a = std::string::npos;
        for (size_t i = 0; i < a->name.length(); ++i) if (isdigit(a->name[i])) { first_digit_a = i; break; }

        size_t first_digit_b = std::string::npos;
        for (size_t i = 0; i < b->name.length(); ++i) if (isdigit(b->name[i])) { first_digit_b = i; break; }

        if (first_digit_a == std::string::npos || first_digit_b == std::string::npos) return a->name < b->name;
        
        std::string prefix_a = a->name.substr(0, first_digit_a);
        std::string prefix_b = b->name.substr(0, first_digit_b);

        if (prefix_a != prefix_b) return prefix_a < prefix_b;
        
        int num_a = std::stoi(a->name.substr(first_digit_a));
        int num_b = std::stoi(b->name.substr(first_digit_b));

        return num_a < num_b; });

        // 7. 遍歷排序後的副本 (sorted_output_nodes) 來進行輸出
        for (const auto &node : sorted_output_nodes)
        {
            const Shape &current_shape = node->module_info->possible_shapes[node->current_shape_idx];
            // 這裡的輸出精度控制邏輯保持不變
            out_file << node->name << " "
                     << fixed << setprecision(4) << node->x << " "
                     << fixed << setprecision(4) << node->y << " ("
                     << fixed << setprecision(2) << current_shape.width << " "
                     << fixed << setprecision(2) << current_shape.height << " "
                     << current_shape.col_multiple << " "
                     << current_shape.row_multiple << ")" << endl;
        }

        out_file.close();
    };
    void packing()
    {
        if (root == nullptr)
            return;

        // 1. 重置所有模組的座標
        for (auto &node : this->nodes)
        {
            node->x = 0.0;
            node->y = 0.0;
        }

        // 2. 每次 packing 前，都重置晶片尺寸和 Contour
        DoublyLinkedList contour;
        current_chip_width = 0.0;
        current_chip_height = 0.0;

        // 3. 從根節點開始進行 DFS Packing
        packing_dfs(root, contour);
    }
    void write_best_tree(const std::string &filename)
    {
        ofstream out_file(filename);
        if (!out_file.is_open())
        {
            cerr << "錯誤: 無法開啟樹結構輸出檔案 " << filename << endl;
            return;
        }

        // 寫入 Graphviz 的標頭
        out_file << "digraph BStarTree {" << endl;
        out_file << "    rankdir=TB;" << endl;
        out_file << "    node [shape=circle, style=filled, fillcolor=lightblue];" << endl;
        out_file << "    edge [fontsize=10];" << endl;

        // 標示出根節點
        if (this->root != nullptr)
        {
            out_file << "    \"" << this->root->name << "\" [fillcolor=salmon];" << endl;
            // 呼叫遞迴函式來產生所有的邊
            generate_edge_list_recursive(this->root, out_file);
        }
        else
        {
            out_file << "    label=\"Tree is empty.\";" << endl;
        }

        // 寫入結尾
        out_file << "}" << endl;

        out_file.close();
    }

private:
    struct Node
    {
        // 模組資訊
        std::string name;
        double width, height;
        int current_shape_idx;
        const Module *module_info; // 新增一個指標，指向原始 Module

        // 座標
        double x, y;

        // B*-Tree 結構指標
        Node *parent;
        Node *left_child;
        Node *right_child;

        // 建構子
        Node(std::string n, double w, double h) : name(n), width(w), height(h), current_shape_idx(0), module_info(nullptr),
                                                  x(0), y(0),
                                                  parent(nullptr), left_child(nullptr), right_child(nullptr) {}
    };

    Node *root;
    std::vector<Node *> nodes;
    Node *best_root;
    std::vector<Node *> best_nodes_state;
    // 新增這兩個成員變數
    double current_chip_width;
    double current_chip_height;

    // --- 實作 Contour 類別 (作為 Placer 的私有內部類別) ---
    struct ContourNode
    {
        // contour 邏輯核心欄位
        string module_id; // 模組 ID (如果需要的話)
        double x_start;   // 線段的起始 x 座標
        double x_end;
        double y; // 線段的高度
    };
    struct ContourDLLNode
    {
        ContourNode data;     // 存放輪廓線段的資料
        ContourDLLNode *prev; // 指向前一個節點
        ContourDLLNode *next; // 指向後一個節點

        // Constructor
        ContourDLLNode(const ContourNode &d) : data(d), prev(nullptr), next(nullptr) {}
    };
    class DoublyLinkedList
    {
    public:
        ContourDLLNode *head;
        ContourDLLNode *tail;

        DoublyLinkedList() : head(nullptr), tail(nullptr) {}

        // Destructor (記得釋放記憶體)
        ~DoublyLinkedList()
        {
            ContourDLLNode *current = head;
            while (current != nullptr)
            {
                ContourDLLNode *next = current->next;
                delete current;
                current = next;
            }
        }

        // 在指定節點之後插入新節點
        void insert_after(ContourDLLNode *node, const ContourNode &data)
        {
            ContourDLLNode *new_node = new ContourDLLNode(data);
            if (node == nullptr)
            { // 插入到開頭
                if (head == nullptr)
                { // 串列為空
                    head = tail = new_node;
                }
                else
                {
                    new_node->next = head;
                    head->prev = new_node;
                    head = new_node;
                }
            }
            else
            { // 在 node 後面插入
                new_node->next = node->next;
                new_node->prev = node;
                if (node->next != nullptr)
                {
                    node->next->prev = new_node;
                }
                else
                {
                    tail = new_node; // 更新 tail
                }
                node->next = new_node;
            }
        }

        // 移除一個節點
        void remove_node(ContourDLLNode *node)
        {
            if (node == nullptr)
                return;

            if (node->prev)
            {
                node->prev->next = node->next;
            }
            else
            { // 移除的是 head
                head = node->next;
            }

            if (node->next)
            {
                node->next->prev = node->prev;
            }
            else
            { // 移除的是 tail
                tail = node->prev;
            }

            delete node;
        }
        // 移除一個節點
        void replace_node(ContourDLLNode *node, const ContourNode &data)
        {
            if (node == nullptr)
            {
                return;
            }
            ContourDLLNode *new_node = new ContourDLLNode(data);

            new_node->next = node->next;
            new_node->prev = node->prev;
            if (node->prev != nullptr)
            {
                node->prev->next = new_node;
            }
            else
            {
                head = new_node; // 更新 head
            }
            if (node->next != nullptr)
            {
                node->next->prev = new_node;
            }
            else
            {
                tail = new_node; // 更新 tail
            }
            delete node;
        }

        bool is_empty() const
        {
            return head == nullptr;
        }
    };
    // --- Packing 演算法的私有遞迴核心 (DFS) ---
    void packing_dfs(Node *node, DoublyLinkedList &contour)
    {
        // Base Case: 如果節點為空，則返回
        if (node == nullptr)
        {
            return;
        }

        // --- 核心邏輯：處理當前節點 (Pre-order Traversal: Node -> Left -> Right) ---

        // 1. 根據父節點決定 X 座標
        if (node->parent == nullptr)
        { // Case: 根節點
            node->x = 0.0;
        }
        else if (node == node->parent->left_child)
        { // Case: 左子節點
            node->x = node->parent->x + node->parent->width;
        }
        else
        { // Case: 右子節點
            node->x = node->parent->x;
        }

        // 2. 根據輪廓線 (Contour) 和父節點決定 Y 座標
        //    (此處假設您已為 std::list<ContourNode> 實現了 get_max_height)
        double y_from_contour = get_max_height_from_contour(contour, node->x, node->x + node->width);

        if (node->parent != nullptr && node == node->parent->right_child)
        {
            // 右子節點必須在父節點的上方
            double y_from_parent = node->parent->y + node->parent->height;
            node->y = std::max(y_from_contour, y_from_parent);
        }
        else
        {
            node->y = y_from_contour;
        }

        // 3. 更新輪廓線和晶片總尺寸
        //    (此處假設您已為 std::list<ContourNode> 實現了 update_contour)
        double new_top_y = node->y + node->height;
        update_contour(contour, node->x, node->x + node->width, new_top_y);

        current_chip_width = std::max(current_chip_width, node->x + node->width);
        current_chip_height = std::max(current_chip_height, node->y + node->height);

        // --- 遞迴呼叫：先走訪左子樹，再走訪右子樹 ---
        packing_dfs(node->left_child, contour);
        packing_dfs(node->right_child, contour);
    }
    double get_max_height_from_contour(DoublyLinkedList &contour, double x_start, double x_end)
    {
        double max_y = 0;
        ContourDLLNode *current = contour.head;

        while (current != nullptr)
        {
            // 檢查是否有重疊
            if (current->data.x_end > x_start && current->data.x_start < x_end)
            {
                // 有重疊，更新最大高度
                max_y = max(max_y, current->data.y);
            }
            current = current->next;
        }
        return max_y;
    }

    // 參數：contour串列, 新模組的x起點, x終點, 和它放置後頂部的y座標
    void update_contour(DoublyLinkedList &contour, double new_x_start, double new_x_end, double new_y)
    {
        // 用於記錄新線段應該插入在哪個節點的後面
        ContourDLLNode *insert_after_node = nullptr;

        ContourDLLNode *current = contour.head;
        while (current != nullptr)
        {
            // 為了安全地刪除節點，預先保存下一個節點的指標
            ContourDLLNode *next_node = current->next;

            double old_x_start = current->data.x_start;
            double old_x_end = current->data.x_end;

            // --- 檢查重疊 ---
            // 如果沒有重疊，且當前節點在新節點的左側，更新 insert_after_node 並繼續
            if (old_x_end <= new_x_start)
            {
                insert_after_node = current;
                current = next_node;
                continue;
            }
            // 如果沒有重疊，且當前節點在新節點的右側，那麼可以提前結束迴圈
            if (old_x_start >= new_x_end)
            {
                break;
            }

            // --- 至此，代表線段有重疊 ---

            // 情況 1: 舊線段的左邊部分沒有被覆蓋，需要保留
            if (old_x_start < new_x_start)
            {
                current->data.x_end = new_x_start; // 縮短舊線段的長度
                insert_after_node = current;       // 新線段應該插在這個線段後面
            }

            // 情況 2: 舊線段的右邊部分沒有被覆蓋，需要保留
            if (old_x_end > new_x_end)
            {
                current->data.x_start = new_x_end; // 縮短舊線段的起點
            }

            // 情況 3: 舊線段被新線段完全覆蓋 (old_x_start >= new_x_start && old_x_end <= new_x_end)
            // 在這種情況下，上面的兩個if都不會執行，所以舊線段的範圍不變，我們需要刪除它。
            // 或者，舊線段的左邊部分被縮短後，右邊又沒超出，那剩下的部分就是被覆蓋的。
            if (old_x_start >= new_x_start && old_x_end <= new_x_end)
            {
                contour.remove_node(current); // 安全地移除節點
            }

            current = next_node;
        }

        // 迴圈結束後，在紀錄的位置插入代表新模組頂部的新線段
        ContourNode new_segment_data = {"", new_x_start, new_x_end, new_y};
        contour.insert_after(insert_after_node, new_segment_data);
    }
    // --- 需要實作 (5): INL 計算 ---
    double calculate_inl()
    {
        if (nodes.empty())
            return 0.0;

        // 步驟 1 & 2: 計算 Bounding Box 和質心
        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::min();
        double y_min = std::numeric_limits<double>::max();
        double y_max = std::numeric_limits<double>::min();

        for (const auto &node : nodes)
        {
            x_min = std::min(x_min, node->x);
            x_max = std::max(x_max, node->x + node->width);
            y_min = std::min(y_min, node->y);
            y_max = std::max(y_max, node->y + node->height);
        }
        double Xc = (x_min + x_max) / 2.0;
        double Yc = (y_min + y_max) / 2.0;

        // 步驟 3: 計算每個 Block 的平方歐幾里得距離
        std::map<std::string, double> squared_distances;
        for (const auto &node : nodes)
        {
            double block_center_x = node->x + node->width / 2.0;
            double block_center_y = node->y + node->height / 2.0;
            squared_distances[node->name] = pow(block_center_x - Xc, 2) + pow(block_center_y - Yc, 2);
        }

        // 步驟 4: 依名稱升序排序
        std::vector<Node *> sorted_nodes = nodes;
        sort(sorted_nodes.begin(), sorted_nodes.end(), [](const Node *a, const Node *b)
             {
        // 找到 a 和 b 名稱中第一個數字的位置
        size_t first_digit_a = std::string::npos;
        for (size_t i = 0; i < a->name.length(); ++i) {
            if (isdigit(a->name[i])) {
                first_digit_a = i;
                break;
            }
        }

        size_t first_digit_b = std::string::npos;
        for (size_t i = 0; i < b->name.length(); ++i) {
            if (isdigit(b->name[i])) {
                first_digit_b = i;
                break;
            }
        }

        // 如果任一名稱中沒有數字，則退回簡單的字串比較
        if (first_digit_a == std::string::npos || first_digit_b == std::string::npos) {
            return a->name < b->name;
        }

        // 提取字母前綴
        std::string prefix_a = a->name.substr(0, first_digit_a);
        std::string prefix_b = b->name.substr(0, first_digit_b);

        // 如果前綴不同，直接比較前綴
        if (prefix_a != prefix_b) {
            return prefix_a < prefix_b;
        }

        // 如果前綴相同，提取數字部分並轉換為整數進行比較
        int num_a = std::stoi(a->name.substr(first_digit_a));
        int num_b = std::stoi(b->name.substr(first_digit_b));

        return num_a < num_b; });
        // 步驟 5: 累積平方距離
        std::vector<double> S_actual;
        double current_sum = 0.0;
        for (const auto &node : sorted_nodes)
        {
            current_sum += squared_distances[node->name];
            S_actual.push_back(current_sum);
        }

        // 步驟 6: 線性迴歸
        double n_sum = 0, s_sum = 0, ns_sum = 0, n2_sum = 0;
        int N = sorted_nodes.size();
        for (int i = 0; i < N; ++i)
        {
            double n = i + 1;
            double s = S_actual[i];
            n_sum += n;
            s_sum += s;
            ns_sum += n * s;
            n2_sum += n * n;
        }

        double a = (N * ns_sum - n_sum * s_sum) / (N * n2_sum - n_sum * n_sum);
        double b = (s_sum - a * n_sum) / N;

        // 步驟 7: 計算 INL
        double inl = 0.0;
        for (int i = 0; i < N; ++i)
        {
            double S_ideal = a * (i + 1) + b;
            inl = std::max(inl, std::abs(S_actual[i] - S_ideal));
        }

        return inl;
    };

    // --- 需要實作 (6): 成本函數 ---
    double calculate_cost()
    {
        double chip_width = this->current_chip_width;   // 直接使用
        double chip_height = this->current_chip_height; // 直接使用

        if (chip_width == 0 || chip_height == 0)
            return std::numeric_limits<double>::max();

        double chip_area = chip_width * chip_height;
        double AR = chip_width / chip_height; // PDF 中的 AR 定義

        // 2. 計算 f(AR)
        double f_AR = 0.0;
        if (AR < 0.5)
        {
            f_AR = 2 * (0.5 - AR);
        }
        else if (AR > 2)
        {
            f_AR = AR - 2;
        }

        // 3. 計算 Cost_AreaAR
        double cost_AreaAR = chip_area * (1.0 + f_AR);

        // 4. 計算 INL
        double inl_val = calculate_inl();

        // 5. 返回加權後的總成本
        // 為了避免因尺度差異導致某項成本被忽略，我們對 inl 做一點縮放
        // 這是一個經驗值，可能需要根據不同 testcase 調整
        return 1 * cost_AreaAR + 10 * inl_val;
    };

    // --- 需要實作 (7): 擾動操作 ---
    void perturb(int type)
    {
        // HINT: 根據傳入的 type，呼叫下面三種操作之一。
        switch (type)
        {
        case 0:
            reshape_node();
            break;
        case 1:
            swap_nodes();
            break;
        case 2:
            delete_and_insert_node();
            break;
        }
    };

    void reshape_node()
    {
        int idx = rand() % nodes.size();
        Node *nd = nodes[idx];
        const auto &shapes = nd->module_info->possible_shapes;
        if (shapes.size() <= 1)
            return;
        int new_idx = rand() % shapes.size();
        nd->current_shape_idx = new_idx;
        // 在 reshape_node() 中
        const auto &new_shape = shapes[new_idx];
        // [修正] 更新總尺寸時，也要乘上倍率
        nd->width = new_shape.width * new_shape.col_multiple;
        nd->height = new_shape.height * new_shape.row_multiple;
    }

    void swap_nodes()
    {
        int n = nodes.size();
        int i = rand() % n, j = rand() % n;
        if (i == j)
            return;
        Node *A = nodes[i], *B = nodes[j];
        // 只交換模組資訊與尺寸，不動樹結構
        std::swap(A->module_info, B->module_info);
        std::swap(A->name, B->name);
        std::swap(A->current_shape_idx, B->current_shape_idx);
        // 重新拿對應 shape 寬高
        // 在 swap_nodes() 中
        {
            const auto &sa = A->module_info->possible_shapes[A->current_shape_idx];
            // [修正] 更新總尺寸時，也要乘上倍率
            A->width = sa.width * sa.col_multiple;
            A->height = sa.height * sa.row_multiple;
        }
        {
            const auto &sb = B->module_info->possible_shapes[B->current_shape_idx];
            // [修正] 更新總尺寸時，也要乘上倍率
            B->width = sb.width * sb.col_multiple;
            B->height = sb.height * sb.row_multiple;
        }
    }
    // 檢查 node_to_check 是否為 ancestor 的子孫
    bool is_descendant(const Node *node_to_check, const Node *ancestor)
    {
        if (ancestor == nullptr || node_to_check == nullptr)
        {
            return false;
        }

        // 使用非遞迴的廣度優先搜尋 (BFS) 或深度優先搜尋 (DFS) 來避免堆疊溢位
        std::queue<const Node *> q;
        if (ancestor->left_child)
            q.push(ancestor->left_child);
        if (ancestor->right_child)
            q.push(ancestor->right_child);

        while (!q.empty())
        {
            const Node *current = q.front();
            q.pop();

            if (current == node_to_check)
            {
                return true; // 找到了，是子孫
            }

            if (current->left_child)
                q.push(current->left_child);
            if (current->right_child)
                q.push(current->right_child);
        }

        return false; // 遍歷完畢，不是子孫
    }
    void delete_and_insert_node()
    {
        if (nodes.size() < 2)
            return;

        // 選一個非 root 節點 A (這部分不變)
        Node *A = nullptr;
        while (true)
        {
            Node *cand = nodes[rand() % nodes.size()];
            if (cand != root)
            {
                A = cand;
                break;
            }
        }

        // 選一個 B，且 B 不是 A 的子孫，且至少有一個空子節點
        std::vector<Node *> cands;
        for (auto nd : nodes)
        {
            // [修改處] 新增 !is_descendant(nd, A) 條件
            if (nd != A && !is_descendant(nd, A) && (nd->left_child == nullptr || nd->right_child == nullptr))
            {
                cands.push_back(nd);
            }
        }

        if (cands.empty())
            return;
        Node *B = cands[rand() % cands.size()];

        // --- 後面的拔除與插入邏輯完全不變 ---
        // 從舊 parent 拔掉 A
        if (A->parent)
        {
            if (A->parent->left_child == A)
                A->parent->left_child = nullptr;
            else
                A->parent->right_child = nullptr;
        }
        A->parent = nullptr;

        // 插到 B 底下
        bool place_left;
        if (B->left_child == nullptr && B->right_child == nullptr)
            place_left = (rand() % 2 == 0);
        else if (B->left_child == nullptr)
            place_left = true;
        else
            place_left = false;

        if (place_left)
            B->left_child = A;
        else
            B->right_child = A;
        A->parent = B;
    }
    // --- 需要實作 (8): 輔助函式 ---
    // --- 修正版: 輔助函式 ---
    // 函式簽名改變，明確傳遞 from_root 和 to_root
    void copy_tree(const std::vector<Node *> &from_nodes, Node *from_root, std::vector<Node *> &to_nodes, Node *&to_root)
    {
        for (Node *node : to_nodes)
        {
            delete node;
        }
        to_nodes.clear();
        to_root = nullptr;

        if (from_nodes.empty())
            return;

        std::map<const Node *, Node *> node_map;

        for (const auto *old_node : from_nodes)
        {
            Node *new_node = new Node(old_node->name, old_node->width, old_node->height);
            new_node->current_shape_idx = old_node->current_shape_idx;
            new_node->module_info = old_node->module_info;
            new_node->x = old_node->x;
            new_node->y = old_node->y;
            to_nodes.push_back(new_node);
            node_map[old_node] = new_node;
        }

        for (const auto *old_node : from_nodes)
        {
            Node *new_node = node_map.at(old_node);
            new_node->parent = old_node->parent ? node_map.at(old_node->parent) : nullptr;
            new_node->left_child = old_node->left_child ? node_map.at(old_node->left_child) : nullptr;
            new_node->right_child = old_node->right_child ? node_map.at(old_node->right_child) : nullptr;
        }

        // 使用傳入的 from_root 來安全地設定新的 to_root
        if (from_root)
        {
            to_root = node_map.at(from_root);
        }
    }

    // 函式簽名改變，明確傳遞 from_root 和 to_root
    void restore_tree(const std::vector<Node *> &from_nodes, Node *from_root, std::vector<Node *> &to_nodes, Node *&to_root)
    {
        // 先釋放 to_nodes (即將被覆蓋的當前狀態) 的記憶體
        for (Node *node : to_nodes)
        {
            delete node;
        }
        to_nodes.clear();
        to_root = nullptr;

        // 重新從 from_nodes 深度複製一份
        // 注意：這裡直接呼叫 copy_tree 會造成重複刪除 to_nodes，所以直接把邏輯寫在這裡
        if (from_nodes.empty())
            return;
        std::map<const Node *, Node *> node_map;
        for (const auto *old_node : from_nodes)
        {
            Node *new_node = new Node(old_node->name, old_node->width, old_node->height);
            new_node->current_shape_idx = old_node->current_shape_idx;
            new_node->module_info = old_node->module_info;
            new_node->x = old_node->x;
            new_node->y = old_node->y;
            to_nodes.push_back(new_node);
            node_map[old_node] = new_node;
        }
        for (const auto *old_node : from_nodes)
        {
            Node *new_node = node_map.at(old_node);
            new_node->parent = old_node->parent ? node_map.at(old_node->parent) : nullptr;
            new_node->left_child = old_node->left_child ? node_map.at(old_node->left_child) : nullptr;
            new_node->right_child = old_node->right_child ? node_map.at(old_node->right_child) : nullptr;
        }
        if (from_root)
        {
            to_root = node_map.at(from_root);
        }
    }
    void generate_edge_list_recursive(const Node *node, ofstream &out)
    {
        if (node == nullptr)
        {
            return;
        }

        // 如果有左子節點，寫入一條邊
        if (node->left_child != nullptr)
        {
            out << "\"" << node->name << "\" -> \"" << node->left_child->name << "\" [label=\"L\"];" << endl;
        }
        // 如果有右子節點，寫入一條邊
        if (node->right_child != nullptr)
        {
            out << "\"" << node->name << "\" -> \"" << node->right_child->name << "\" [label=\"R\"];" << endl;
        }

        // 遞迴處理
        generate_edge_list_recursive(node->left_child, out);
        generate_edge_list_recursive(node->right_child, out);
    }
};

std::vector<Module> parse_block_file(const std::string &filename)
{
    std::vector<Module> modules;
    std::ifstream input_file(filename);

    if (!input_file.is_open())
    {
        std::cerr << "Error: 無法開啟檔案 " << filename << std::endl;
        return modules;
    }

    std::string line;
    while (std::getline(input_file, line))
    {
        if (line.empty())
        {
            continue;
        }

        std::stringstream line_stream(line);
        Module current_module;

        line_stream >> current_module.name;
        char parenthesis;
        while (line_stream >> parenthesis && parenthesis == '(')
        {
            Shape shape;
            line_stream >> shape.width >> shape.height >> shape.col_multiple >> shape.row_multiple;
            line_stream >> parenthesis;
            current_module.possible_shapes.push_back(shape);
        }
        modules.push_back(current_module);
    }

    input_file.close();
    return modules;
}
int main(int argc, char *argv[])
{
    srand(time(0));
    if (argc < 3)
    {
        std::cerr << "用法: " << argv[0] << " <input_file.block> <output_file.output>" << std::endl;
        return 1;
    }

    std::string input_filename = argv[1];
    std::string output_filename = argv[2];

    std::vector<Module> my_modules = parse_block_file(input_filename);

    std::cout << "成功解析 " << my_modules.size() << " 個模組。" << std::endl;

    Placer my_placer(my_modules);
    my_placer.solve();
    my_placer.write_output(output_filename);

    // [新增] 呼叫新的函式來輸出樹結構
    my_placer.write_best_tree("best_tree_structure.txt");

    std::cout << "Packing 完成，請檢查輸出檔案 " << output_filename << std::endl;
    // [新增] 提示使用者有額外的輸出檔案
    std::cout << "最佳 B*-Tree 結構已輸出至 best_tree_structure.txt" << std::endl;

    return 0;
}
// int main(int argc, char *argv[])
// {
//     srand(time(0));
//     if (argc < 3)
//     {
//         std::cerr << "用法: " << argv[0] << " <input_file.block> <output_file.output>" << std::endl;
//         return 1;
//     }

//     std::string input_filename = argv[1];
//     std::string output_filename = argv[2];

//     std::vector<Module> my_modules = parse_block_file(input_filename);

//     std::cout << "成功解析 " << my_modules.size() << " 個模組。" << std::endl;

//     Placer my_placer(my_modules);
//     my_placer.solve();
//     my_placer.write_output(output_filename);

//     std::cout << "Packing 完成，請檢查輸出檔案 " << output_filename << std::endl;

//     return 0;
// }