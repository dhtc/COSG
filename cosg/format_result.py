def format_degs(adata, key, cal_logFC = True, wide = True):
	result = adata.uns[key]
	groups = result['names'].dtype.names
	
	if wide:
		if cal_logFC:
			out_n = ['gene', 'logFC', 'score']
			degs = pd.DataFrame(
				{group + '_' + out_n[i]: result[key][group]
				for group in groups for i, key in enumerate(['names', 'logfoldchanges', 'scores'])})
		else:
			out_n = ['gene', 'score']
			degs = pd.DataFrame(
				{group + '_' + out_n[i]: result[key][group]
				for group in groups for i, key in enumerate(['names', 'scores'])})
	else:
		if cal_logFC:
			degs = pd.concat([pd.DataFrame((result['names'][group], result['logfoldchanges'][group], result['scores'][group], [group] * (result['scores'][group].shape[0]))).T for group in groups], axis=0)
			degs.columns = ['gene', 'logFC', 'score', 'group']
		else:
			degs = pd.concat([pd.DataFrame((result['names'][group], result['scores'][group], [group] * (result['scores'][group].shape[0]))).T for group in groups], axis=0)
			degs.columns = ['gene', 'score', 'group']
	return degs
	
