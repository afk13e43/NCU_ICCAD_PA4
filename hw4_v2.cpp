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
#include <random>
#include <ctime> // 確保 time() 被宣告

using namespace std;

// 用於儲存單一形態的幾何資訊
struct Shape
{
    double width;
    double height;
    int col_multiple;
    int row_multiple;
};

// 用於儲存一個裝置 (模組) 的資訊
struct Module
{
    std::string name;
    std::vector<Shape> possible_shapes;
};

class Placer
{
public:
    // B*-Tree 節點定義
    struct Node
    {
        // 模組資訊
        std::string name;
        double width, height;
        int current_shape_idx;
        const Module *module_info;

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

        // 複製建構子，用於深度複製狀態
        Node(const Node &other) : name(other.name), width(other.width), height(other.height),
                                  current_shape_idx(other.current_shape_idx), module_info(other.module_info),
                                  x(other.x), y(other.y),
                                  parent(nullptr), left_child(nullptr), right_child(nullptr) {}
    };

    Placer(const std::vector<Module> &modules)
    {
        root = nullptr;

        for (const auto &module : modules)
        {
            const Shape &initial_shape = module.possible_shapes[0];
            Node *new_node = new Node(module.name, initial_shape.width, initial_shape.height);
            new_node->module_info = &module;
            this->nodes.push_back(new_node);
        }

        if (nodes.empty())
            return;

        root = nodes[0];

        for (size_t i = 1; i < nodes.size(); ++i)
        {
            Node *new_node = nodes[i];
            bool placed = false;
            while (!placed)
            {
                int parent_idx = rand() % i;
                Node *parent_node = nodes[parent_idx];

                if (rand() % 2 == 0)
                {
                    if (parent_node->left_child == nullptr)
                    {
                        parent_node->left_child = new_node;
                        new_node->parent = parent_node;
                        placed = true;
                    }
                }
                else
                {
                    if (parent_node->right_child == nullptr)
                    {
                        parent_node->right_child = new_node;
                        new_node->parent = parent_node;
                        placed = true;
                    }
                }
            }
        }
    }

    // --- 新增: 解構子來釋放記憶體 ---
    ~Placer()
    {
        for (auto node : nodes)
        {
            delete node;
        }
        for (auto node : best_nodes_state)
        {
            delete node;
        }
    }

    void solve()
    {
        // 1. 初始參數設定
        double T = 10000.0;              // 初始溫度
        double T_min = 0.1;              // 終止溫度
        double cooling_rate = 0.995;     // 降溫速率 (稍微減慢，效果較好)
        int max_iter = nodes.size() * 5; // 每個溫度下的迭代次數，與問題規模相關

        // 2. 產生初始解並計算成本
        packing();
        double current_cost = calculate_cost();
        double best_cost = current_cost;
        copy_tree(nodes, best_nodes_state); // 備份初始最佳解
        cout << "Initial Cost: " << best_cost << endl;

        // 3. 模擬退火主迴圈
        while (T > T_min)
        {
            for (int i = 0; i < max_iter; ++i)
            {
                // a. 儲存當前狀態 (以便退回)
                std::vector<Node *> backup_nodes_state;
                copy_tree(nodes, backup_nodes_state);
                Node *backup_root = root;

                // b. 擾動 (Perturb)，產生新解
                int move_type = rand() % 3; // 修正: 擾動類型是 0, 1, 2
                perturb(move_type);

                // c. 計算新解的成本
                packing();
                double new_cost = calculate_cost();

                // d. 判斷是否接受新解
                double delta_cost = new_cost - current_cost;
                if (delta_cost < 0)
                {
                    current_cost = new_cost;
                    if (new_cost < best_cost)
                    {
                        best_cost = new_cost;
                        copy_tree(nodes, best_nodes_state); // 更新全域最佳解
                        cout << "New Best Cost: " << best_cost << " at T=" << T << endl;
                    }
                    // 清理備份，因為我們接受了新狀態
                    for (auto n : backup_nodes_state)
                        delete n;
                }
                else
                {
                    double accept_prob = exp(-delta_cost / T);
                    if ((double)rand() / RAND_MAX < accept_prob)
                    {
                        current_cost = new_cost;
                        // 清理備份
                        for (auto n : backup_nodes_state)
                            delete n;
                    }
                    else
                    {
                        // 不接受，還原到擾動前的狀態
                        restore_tree(backup_nodes_state, nodes);
                        root = this->nodes[0]; // 確保 root 指標正確
                    }
                }
            }
            // e. 降溫
            T *= cooling_rate;
        }

        // 4. 迴圈結束，還原到找到的最佳解
        restore_tree(best_nodes_state, nodes);
        root = this->nodes[0];
        packing();
        cout << "SA finished. Final Best Cost: " << best_cost << endl;
    }

    void write_output(const std::string &filename)
    {
        packing(); // 確保輸出的是最後狀態的 packing 結果

        double chip_width = this->current_chip_width;
        double chip_height = this->current_chip_height;
        double chip_area = chip_width * chip_height;
        double inl = calculate_inl();

        ofstream out_file(filename);
        if (!out_file.is_open())
        {
            cerr << "錯誤: 無法開啟輸出檔案 " << filename << endl;
            return;
        }

        out_file << fixed;
        out_file << setprecision(4) << chip_area << endl;
        out_file << setprecision(2) << chip_width << " " << chip_height << endl;
        out_file << setprecision(2) << inl << endl;

        // 排序節點以符合某些驗證器的要求 (可選)
        sort(nodes.begin(), nodes.end(), [](const Node *a, const Node *b)
             { return a->name < b->name; });

        for (const auto &node : nodes)
        {
            const Shape &current_shape = node->module_info->possible_shapes[node->current_shape_idx];
            out_file << node->name << " "
                     << node->x << " "
                     << node->y << " ("
                     << current_shape.width << " "
                     << current_shape.height << " "
                     << current_shape.col_multiple << " "
                     << current_shape.row_multiple << ")" << endl;
        }
        out_file.close();
    };
    void packing()
    {
        if (root == nullptr)
            return;

        // 每次 packing 前重置狀態
        Contour contour;
        current_chip_width = 0.0;
        current_chip_height = 0.0;

        // 呼叫遞迴函式開始 packing
        packing_recursive(root, contour);
    }

private:
    Node *root;
    std::vector<Node *> nodes;
    std::vector<Node *> best_nodes_state;
    double current_chip_width;
    double current_chip_height;

    // Contour 類別 (與你原本的實作相同，此處省略以節省空間)
    // ... Contour class ...
    struct ContourNode
    {
        double x;
        double y;
        ContourNode *prev;
        ContourNode *next;
    };
    // Placer class 內部...

    class Contour
    {
    public:
        // map<x_coordinate, y_height>
        // 代表從 x 座標開始，高度變為 y
        std::map<double, double> segments;

        Contour()
        {
            // 初始輪廓：從 x=0 開始，高度為 0
            segments[0.0] = 0.0;
        }

        double get_max_height(double x_start, double x_end)
        {
            double max_y = 0.0;

            // upper_bound 找到第一個 > x_start 的點
            // 我們需要從它前一個點開始看起
            auto it = segments.upper_bound(x_start);
            if (it != segments.begin())
            {
                it--;
            }

            while (it != segments.end() && it->first < x_end)
            {
                max_y = std::max(max_y, it->second);
                it++;
            }
            return max_y;
        }

        void update(double x, double y, double width, double height)
        {
            double x_end = x + width;
            double new_top_y = y + height;

            // 找到 x_end 位置原本的高度，以便稍後恢復
            double original_y_at_x_end = get_max_height(x_end, x_end + 1e-9);

            // 刪除 [x, x_end) 區間內所有舊的輪廓點
            auto it = segments.upper_bound(x);
            while (it != segments.end() && it->first < x_end)
            {
                it = segments.erase(it);
            }

            // 插入新的、更高的輪廓點
            segments[x] = new_top_y;

            // 在 x_end 處恢復輪廓線
            // 只有當 x_end 不在 map 中，或在 map 中但高度需要更新時，才插入
            if (segments.find(x_end) == segments.end() || segments[x_end] < original_y_at_x_end)
            {
                segments[x_end] = original_y_at_x_end;
            }

            // 清理冗餘的點：如果一個點的高度和前一個點相同，則可以移除
            it = segments.find(x);
            if (it != segments.begin() && it->second == std::prev(it)->second)
            {
                segments.erase(it);
            }
            it = segments.find(x_end);
            if (it != segments.begin() && it->second == std::prev(it)->second)
            {
                segments.erase(it);
            }
        }
    };

    void packing_recursive(Node *node, Contour &contour)
    {
        if (node == nullptr)
            return;

        // 根據父節點和輪廓線，決定目前節點的 x, y
        if (node == root)
        {
            node->x = 0;
            node->y = 0;
        }
        else
        {
            Node *parent = node->parent;
            if (parent->left_child == node)
            { // 是左子節點
                node->x = parent->x + parent->width;
                node->y = contour.get_max_height(node->x, node->x + node->width);
            }
            else
            { // 是右子節點
                node->x = parent->x;
                double y_from_contour = contour.get_max_height(node->x, node->x + node->width);
                double y_from_parent = parent->y + parent->height; // B*Tree 的約束
                node->y = std::max(y_from_contour, y_from_parent);
            }
        }

        // 更新輪廓線和晶片總尺寸
        contour.update(node->x, node->y, node->width, node->height);
        current_chip_width = std::max(current_chip_width, node->x + node->width);
        current_chip_height = std::max(current_chip_height, node->y + node->height);

        // 遞迴處理子節點 (Pre-order Traversal)
        packing_recursive(node->left_child, contour);
        packing_recursive(node->right_child, contour);
    }
    double calculate_inl()
    {
        // 你的 INL 實作是正確的，直接使用
        if (nodes.empty())
            return 0.0;

        double Xc = current_chip_width / 2.0;
        double Yc = current_chip_height / 2.0;

        std::map<std::string, double> squared_distances;
        for (const auto &node : nodes)
        {
            double block_center_x = node->x + node->width / 2.0;
            double block_center_y = node->y + node->height / 2.0;
            squared_distances[node->name] = pow(block_center_x - Xc, 2) + pow(block_center_y - Yc, 2);
        }

        std::vector<Node *> sorted_nodes = nodes;
        sort(sorted_nodes.begin(), sorted_nodes.end(), [](const Node *a, const Node *b)
             { return a->name < b->name; });

        std::vector<double> S_actual;
        double current_sum = 0.0;
        for (const auto &node : sorted_nodes)
        {
            current_sum += squared_distances[node->name];
            S_actual.push_back(current_sum);
        }

        if (S_actual.size() <= 1)
            return 0.0;

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

        double den = (N * n2_sum - n_sum * n_sum);
        if (abs(den) < 1e-9)
            return 0.0; // 避免除以零

        double a = (N * ns_sum - n_sum * s_sum) / den;
        double b = (s_sum - a * n_sum) / N;

        double inl = 0.0;
        for (int i = 0; i < N; ++i)
        {
            double S_ideal = a * (i + 1) + b;
            inl = std::max(inl, std::abs(S_actual[i] - S_ideal));
        }

        return inl;
    };

    double calculate_cost()
    {
        double chip_width = this->current_chip_width;
        double chip_height = this->current_chip_height;

        if (chip_width < 1e-6 || chip_height < 1e-6)
            return std::numeric_limits<double>::max();

        double chip_area = chip_width * chip_height;
        double AR = max(chip_width / chip_height, chip_height / chip_width);

        double f_AR = 0.0;
        if (AR > 2.0)
            f_AR = AR - 2.0;

        double cost_AreaAR = chip_area * (1.0 + f_AR);
        double inl_val = calculate_inl();

        // --- 修正: 調整權重 ---
        // Area 和 INL 的數值級別差異可能很大，需要用權重平衡
        // alpha 通常為 1，beta 需要實驗調整，這裡給一個經驗值
        double alpha = 1.0;
        double beta = 1e-7;

        return alpha * cost_AreaAR + beta * inl_val;
    };

    void perturb(int type)
    {
        switch (type)
        {
        case 0:
            reshape_node();
            break;
        case 1:
            swap_nodes_content();
            break;
        case 2:
            delete_and_insert_leaf_node();
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
        nd->width = shapes[new_idx].width;
        nd->height = shapes[new_idx].height;
    }

    void swap_nodes_content()
    {
        if (nodes.size() < 2)
            return;
        int i = rand() % nodes.size();
        int j = rand() % nodes.size();
        if (i == j)
            return;
        Node *A = nodes[i], *B = nodes[j];

        // 只交換內容，不改變樹的拓撲結構，這是最簡單安全的交換
        std::swap(A->name, B->name);
        std::swap(A->module_info, B->module_info);
        std::swap(A->current_shape_idx, B->current_shape_idx);
        std::swap(A->width, B->width);
        std::swap(A->height, B->height);
    }

    // --- 修正: delete_and_insert_node ---
    // 原本的實作會造成孤兒節點。修正為只移動葉節點，這是最簡單且保證樹結構完整的作法。
    void delete_and_insert_leaf_node()
    {
        if (nodes.size() < 2)
            return;

        // 1. 找到一個非 root 的葉節點 (leaf node) 來移動
        std::vector<Node *> movable_nodes;
        for (auto node : nodes)
        {
            if (node != root && node->left_child == nullptr && node->right_child == nullptr)
            {
                movable_nodes.push_back(node);
            }
        }
        if (movable_nodes.empty())
            return; // 沒有可移動的葉節點
        Node *node_to_move = movable_nodes[rand() % movable_nodes.size()];

        // 2. 從樹中安全地移除它
        Node *parent = node_to_move->parent;
        if (parent->left_child == node_to_move)
        {
            parent->left_child = nullptr;
        }
        else
        {
            parent->right_child = nullptr;
        }
        node_to_move->parent = nullptr;

        // 3. 找到一個可以安插的新位置 (目標節點)
        std::vector<Node *> target_nodes;
        for (auto node : nodes)
        {
            // 不能是自己，且至少有一個空位
            if (node != node_to_move && (node->left_child == nullptr || node->right_child == nullptr))
            {
                target_nodes.push_back(node);
            }
        }
        if (target_nodes.empty())
        { // 如果找不到，就把它插回原位 (操作失敗)
            if (parent->left_child == nullptr)
                parent->left_child = node_to_move;
            else
                parent->right_child = node_to_move;
            node_to_move->parent = parent;
            return;
        }
        Node *target_node = target_nodes[rand() % target_nodes.size()];

        // 4. 安插到目標節點下
        if (target_node->left_child == nullptr && target_node->right_child == nullptr)
        {
            if (rand() % 2 == 0)
                target_node->left_child = node_to_move;
            else
                target_node->right_child = node_to_move;
        }
        else if (target_node->left_child == nullptr)
        {
            target_node->left_child = node_to_move;
        }
        else
        {
            target_node->right_child = node_to_move;
        }
        node_to_move->parent = target_node;
    }

    // --- 已實作: copy_tree ---
    void copy_tree(const std::vector<Node *> &source_nodes, std::vector<Node *> &dest_nodes)
    {
        // 清理目標 vector，避免記憶體洩漏
        for (auto node : dest_nodes)
        {
            delete node;
        }
        dest_nodes.clear();

        if (source_nodes.empty())
            return;

        map<const Node *, Node *> node_map; // <原始節點, 複製節點>

        // 1. 深度複製所有節點，但不設定指標
        for (const auto &src_node : source_nodes)
        {
            Node *new_node = new Node(*src_node); // 使用複製建構子
            dest_nodes.push_back(new_node);
            node_map[src_node] = new_node;
        }

        // 2. 根據 map 重建指標連結
        for (const auto &src_node : source_nodes)
        {
            Node *new_node = node_map[src_node];
            if (src_node->parent)
            {
                new_node->parent = node_map[src_node->parent];
            }
            if (src_node->left_child)
            {
                new_node->left_child = node_map[src_node->left_child];
            }
            if (src_node->right_child)
            {
                new_node->right_child = node_map[src_node->right_child];
            }
        }
    }

    // --- 已實作: restore_tree ---
    void restore_tree(const std::vector<Node *> &backup_nodes, std::vector<Node *> &active_nodes)
    {
        // 這個函式是用備份的狀態，覆蓋掉當前的活動狀態
        // backup_nodes 是 SA 迴圈中建立的臨時備份，用完需要刪除
        // active_nodes 是 Placer class 的 this->nodes

        map<string, Node *> active_node_map;
        for (auto node : active_nodes)
        {
            active_node_map[node->name] = node;
        }

        map<string, const Node *> backup_node_map;
        for (const auto node : backup_nodes)
        {
            backup_node_map[node->name] = node;
        }

        // 1. 覆蓋內容
        for (const auto &backup_node : backup_nodes)
        {
            Node *active_node = active_node_map[backup_node->name];
            active_node->width = backup_node->width;
            active_node->height = backup_node->height;
            active_node->current_shape_idx = backup_node->current_shape_idx;
        }

        // 2. 覆蓋指標
        for (const auto &backup_node : backup_nodes)
        {
            Node *active_node = active_node_map[backup_node->name];

            if (backup_node->parent)
            {
                active_node->parent = active_node_map[backup_node->parent->name];
            }
            else
            {
                active_node->parent = nullptr;
                this->root = active_node; // 更新 Placer 的 root 指標
            }

            if (backup_node->left_child)
            {
                active_node->left_child = active_node_map[backup_node->left_child->name];
            }
            else
            {
                active_node->left_child = nullptr;
            }

            if (backup_node->right_child)
            {
                active_node->right_child = active_node_map[backup_node->right_child->name];
            }
            else
            {
                active_node->right_child = nullptr;
            }
        }

        // 3. 清理臨時備份的記憶體
        for (auto n : backup_nodes)
            delete n;
    }
};

std::vector<Module> parse_block_file(const std::string &filename)
{
    // 你的 parser 實作是正確的
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
            continue;

        std::stringstream line_stream(line);
        Module current_module;
        line_stream >> current_module.name;

        char parenthesis;
        while (line_stream >> parenthesis && parenthesis == '(')
        {
            Shape shape;
            line_stream >> shape.width >> shape.height >> shape.col_multiple >> shape.row_multiple;
            line_stream >> parenthesis; // ')'
            current_module.possible_shapes.push_back(shape);
        }
        modules.push_back(current_module);
    }

    input_file.close();
    return modules;
}

int main(int argc, char *argv[])
{
    // --- 修正: main 函式參數解析 ---
    if (argc != 5)
    {
        cerr << "用法: " << argv[0] << " -i <input_file.block> -o <output_file.output>" << endl;
        return 1;
    }
    string input_filename, output_filename;
    for (int i = 1; i < argc; ++i)
    {
        string arg = argv[i];
        if (arg == "-i" && i + 1 < argc)
        {
            input_filename = argv[++i];
        }
        else if (arg == "-o" && i + 1 < argc)
        {
            output_filename = argv[++i];
        }
    }
    if (input_filename.empty() || output_filename.empty())
    {
        cerr << "錯誤: 未指定輸入或輸出檔案。" << endl;
        return 1;
    }

    srand(time(0));

    std::vector<Module> my_modules = parse_block_file(input_filename);
    if (my_modules.empty())
    {
        return 1;
    }
    std::cout << "成功解析 " << my_modules.size() << " 個模組。" << std::endl;

    Placer my_placer(my_modules);
    my_placer.solve();
    my_placer.write_output(output_filename);

    std::cout << "Packing 完成，請檢查輸出檔案 " << output_filename << std::endl;

    return 0;
}