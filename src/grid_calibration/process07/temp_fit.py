from .__init__ import tb_op

data_all = tb_op.load_data()
tb_op.temp_bias_fit(data_all)