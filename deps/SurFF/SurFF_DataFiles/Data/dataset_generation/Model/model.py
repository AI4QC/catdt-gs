from typing import Optional

import torch
import torch.nn.functional as F
import torch.nn as nn
from torch.nn import init
from torch import Tensor
from torch_geometric.nn import GCNConv, global_mean_pool, \
    GATConv, GATv2Conv, global_add_pool, TopKPooling, SAGPooling, MessagePassing

from torch_geometric.utils import softmax
from torch_geometric.typing import OptTensor


class GraphMHSA(MessagePassing):
    def __init__(self, in_channels, out_channels, num_heads, edge_dim, dropout_prob=0):
        super(GraphMHSA, self).__init__(node_dim=0, aggr='add')
        self.num_heads = num_heads
        self.head_dim = out_channels // num_heads
        self.scale = self.head_dim ** 0.5
        self.out_channels = out_channels

        self.W_Q = torch.nn.Linear(in_channels, out_channels)
        self.W_K = torch.nn.Linear(in_channels, out_channels)
        self.W_V = torch.nn.Linear(in_channels, out_channels)
        self.W_O = nn.Linear(out_channels, out_channels)

        self.dropout_prob = dropout_prob

    def forward(self, x, edge_index, edge_attr):
        out = self.propagate(edge_index, x=x, edge_attr=edge_attr)
        return out

    def message(self, x_i: Tensor, x_j: Tensor, edge_attr: Tensor, index: Tensor, ptr: OptTensor,
                size_i: Optional[int]) -> Tensor:
        # edge_attr = self.W_E(edge_attr)  # [edge, head]
        # edge_attr = edge_attr.repeat(1, self.num_heads)  # [edge, head]

        Q_i = self.W_Q(x_i).view(-1, self.num_heads, self.head_dim)  # [edge, head, dim]
        K_j = self.W_K(x_j).view(-1, self.num_heads, self.head_dim)
        V_j = self.W_V(x_j).view(-1, self.num_heads, self.head_dim)

        attention_scores = torch.sum(Q_i * K_j, dim=-1)
        attention_scores = attention_scores / self.scale + edge_attr  # [edge, head]

        attention_scores = softmax(attention_scores, index, ptr, size_i)
        # attention_scores = F.dropout(attention_scores, self.dropout_prob, self.training)

        messages = V_j * attention_scores.unsqueeze(-1)

        return messages

    def update(self, attention_output):
        # Concatenate and apply final linear projection

        attention_output = self.concat_heads(attention_output)
        return self.W_O(attention_output)

    def concat_heads(self, x):
        return x.view(-1, self.out_channels)


class Graphconv(nn.Module):
    def __init__(self, conv_feature, conv_channels, edge_dim, dropout=0):
        super().__init__()
        self.conv_layers = GraphMHSA(conv_feature, conv_feature, conv_channels, edge_dim)

        self.ln1 = nn.LayerNorm(conv_feature)
        self.ln2 = nn.LayerNorm(conv_feature)

        self.ffn1 = nn.Linear(conv_feature, conv_feature)
        self.ffn2 = nn.Linear(conv_feature, conv_feature)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, edge_index, edge_attr):
        res = x
        x = self.ln1(x)
        x = self.conv_layers(x=x, edge_index=edge_index, edge_attr=edge_attr)  # Add the conv output to x with residual
        x = self.dropout(x)  # Apply dropout
        x = res + x

        res = x
        x = self.ln2(x)
        x = F.gelu(self.ffn1(x))
        x = self.dropout(x)
        x = self.ffn2(x)
        x = self.dropout(x)
        x = res + x
        return x


@torch.jit.script
def gaussian(x, mean, std):
    pi = 3.14159
    a = (2 * pi) ** 0.5
    return torch.exp(-0.5 * (((x - mean) / std) ** 2)) / (a * std)


class GaussianLayer(nn.Module):
    def __init__(self, K=128, edge_types=1024):
        super().__init__()
        self.K = K
        self.means = nn.Embedding(1, K)
        self.stds = nn.Embedding(1, K)
        self.mul = nn.Embedding(edge_types, 1)
        self.bias = nn.Embedding(edge_types, 1)
        # nn.init.uniform_(self.means.weight, 0, 6)
        # nn.init.uniform_(self.stds.weight, 0, 0.5)

        self.means.weight = nn.Parameter(torch.linspace(0, 6, K).view(1, -1), requires_grad=False)
        self.stds.weight = nn.Parameter((torch.zeros(K) + 0.1).view(1, -1), requires_grad=False)

        nn.init.constant_(self.bias.weight, 0)
        nn.init.constant_(self.mul.weight, 0)

    def forward(self, x, edge_types):
        # mul = self.mul(edge_types) + 1.0
        # bias = self.bias(edge_types)
        # x = mul * x + bias
        x = x.expand(-1, self.K)
        mean = self.means.weight.float().view(-1)
        std = self.stds.weight.float().view(-1).abs() + 1e-5

        return gaussian(x.float(), mean, std).type_as(self.mul.weight)


class GCN(torch.nn.Module):
    def __init__(self, atom_feature=92, conv_feature=512, conv_block=4, fc_feature=256, out_feature=1, dropout=0,
                 num_conv_layers=12, num_final_conv_layers=0, num_fc_layers=2, conv_channels=4, edge_dim=61,
                 layer_fea=False):
        super().__init__()
        # Define the dropout layer with probability dropout
        self.dropout = nn.Dropout(p=dropout)

        # Define a linear embedding layer to convert atom feature to conv feature
        self.embed = nn.Embedding(atom_feature, conv_feature)
        self.cent_enc_proj = nn.Linear(edge_dim, conv_feature)
        # self.edge_embed = nn.Linear(1, edge_dim)
        # Define num_conv_layers number of GATv2Conv layers as a ModuleList

        assert (conv_feature % conv_channels) == 0, AssertionError(
            "num_conv_feature must be divisible by num_conv_channels")

        # self.conv_layers = nn.ModuleList([GATv2Conv(conv_feature, int(conv_feature / conv_channels),
        #                                             edge_dim=edge_dim, heads=conv_channels,
        #                                             dropout=dropout) for _ in range(num_conv_layers)])

        self.edge_embed = GaussianLayer(K=edge_dim, edge_types=100 * 101)
        self.W_E = nn.Linear(edge_dim, conv_channels, bias=False)

        self.conv_layers = nn.ModuleList([Graphconv(conv_feature, conv_channels, edge_dim, dropout=dropout)
                                          for _ in range(num_conv_layers)])

        self.conv_block = conv_block

        ## Define the final conv layers

        self.conv_final = nn.ModuleList([Graphconv(conv_feature, conv_channels, edge_dim, dropout=dropout)
                                         for _ in range(num_final_conv_layers)])

        self.final_ln = nn.LayerNorm(conv_feature)

        self.fc_layers = nn.ModuleList([nn.Linear(fc_feature, fc_feature)] * num_fc_layers)

        # Define the output layer with fc_feature input and out_feature output
        self.fc_out = nn.Linear(fc_feature, out_feature)

    def forward(self, data):
        # Unpack the input data
        x, edge_attr, edge_index, batch, tags, area = data.x, data.edge_attr, data.edge_index, data.batch, data.tags, data.area
        x = x.to(torch.long)  # [N]

        # calculate edge_feature

        edge_attr = self.edge_embed(edge_attr, x[edge_index[0]] * 100 + x[edge_index[1]])  # [N, gaussian_kernel]

        # get centrality encoding of atoms
        # cent_enc = []
        # for i in range(len(x)):
        #     edge_index_i = torch.where(edge_index[0] == i)
        #     encoding_i = torch.sum(edge_attr[edge_index_i], dim=0)
        #     cent_enc.append(encoding_i)
        # cent_enc = torch.vstack(cent_enc) #[N, gaussian_kernel]

        # cent_enc = torch.zeros((len(x), edge_attr.size(1)), device=edge_attr.device)
        # cent_enc.index_add_(0, edge_index[0], edge_attr)

        # cent_enc = self.cent_enc_proj(cent_enc) # [N, dim]

        edge_attr = self.W_E(edge_attr)  # [edge, heads]

        # Apply linear embedding to x to convert it to the conv_feature size
        fea = self.embed(x)
        fea = self.dropout(fea)

        # gcn_fea.append(x)

        # Apply num_conv_layers number of GATv2Conv layers to x
        for _ in range(self.conv_block):
            for i, conv in enumerate(self.conv_layers):
                fea = conv(fea, edge_index, edge_attr)

        for i, conv in enumerate(self.conv_final):
            fea = conv(fea, edge_index, edge_attr)

        fea = self.final_ln(fea)

        # Apply num_fc_layers number of fully connected layers to x
        for fc in self.fc_layers:
            fea = fc(fea)  # Apply fully connected layer
            fea = F.gelu(fea)  # Apply ReLU activation
            fea = self.dropout(fea)  # Apply dropout

        # Apply the output layer to x to get the final prediction

        fea = self.fc_out(fea)

        fea = fea * (tags > 0).unsqueeze(-1).float()

        # num = [torch.sum(tags[data.ptr[i]:data.ptr[i + 1]]) for i in range(len(data.ptr) - 1)]
        # num = torch.stack(num).unsqueeze(-1)

        fea = global_add_pool(fea, batch)
        out = fea

        return out
