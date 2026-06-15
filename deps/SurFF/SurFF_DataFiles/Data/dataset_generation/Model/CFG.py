class CFG:
    """default configurations of model"""

    def __init__(self):
        self.trainset = r""
        self.testset = r""
        self.valset = None

        self.device = "cuda"
        self.batch_size = 64
        self.accumulation_batch = 4

        self.atom_feature = 92
        self.num_edge_fea = 128
        self.conv_feature = 1024
        self.conv_channels = 8
        self.fc_feature = 64
        self.out_feature = 1
        self.num_conv_layers = 3
        self.num_fc_layers = 2
        self.conv_block = 1
        self.layer_fea = True

        self.epoch = 200

        self.lr = 1e-4
        self.weight_decay = 1e-4
        self.lr_patience = 3
        self.gamma = 0.5
        self.dropout = 0

        self.print_freq = 1000
        self.test_freq = 10000
        self.save_freq = 10000

        self.step = 0
        self.save_dir = "test"  # this will create a new folder under the same folder as Run.py
        self.save_name = "test.pth"

        self.model_path = None # this will load an existing model
