import pickle

wulff_pth = r"results/element_exp2&3/wulff_results.pkl"

with open(wulff_pth, 'rb') as f:
    results = pickle.load(f)


