class CustomError(BaseException):

    def __init__(self,msg=None):
        super().__init__(msg)

def validate_input_data(df):
    expected_col_names = ['protein','motif','log2FC','pval']
    column_names = df.columns.to_list()
    if 'protein' not in column_names:
        raise CustomError(msg="Input file error:: Column, 'protein' is missing")
    elif 'motif' not in column_names:
        raise CustomError(msg="Input file error:: Column, 'motif' is missing")
    if len(column_names) > 2:
        for col_name in column_names:
            if col_name not in expected_col_names:
                 raise CustomError(msg=f"Input file error:: Incorrect column {col_name}")
        