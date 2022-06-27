samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["unit"], dtype=str)