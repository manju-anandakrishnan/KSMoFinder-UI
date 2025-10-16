from st_aggrid import GridOptionsBuilder, AgGrid
from st_link_analysis import NodeStyle, EdgeStyle

def build_grid(df,key):
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_selection(selection_mode="single", use_checkbox=True)
    gb.configure_pagination(enabled=True)

    return AgGrid(df,                  
                height=min(27+len(df)*35, 400),                
                gridOptions=gb.build(),
                enable_enterprise_modules=False,              
                fit_columns_on_grid_load=True,
                key=str(key),
                custom_css={
                    "#gridToolBar": {
                        "padding-bottom": "0px !important",
                    }
                })


def get_graph_elements(triples):
    elements = {'nodes':[], 'edges':[]}
    for triple in triples:
        head, rel, tail = triple.get()
        head_node = {'data':{'id':head.get_id(),'name':head.get_name(),'label':head.get_category()}}
        tail_node = {'data':{'id':tail.get_id(),'name':tail.get_name(),'label':tail.get_category()}}
        elements['nodes'].append(head_node)
        elements['nodes'].append(tail_node)

        ht_edge = {'data':{'id':rel.id,'label':rel.get_label(),'source':head.get_id(),'target':tail.get_id()}}
        elements['edges'].append(ht_edge)
    return elements

def get_edge_labels(edges):
    edge_labels = set()
    for edge in edges:
        edge_labels.add(edge['data']['label'])
    return edge_labels

def get_node_labels(nodes):
    node_categories = set()
    for node in nodes:
        node_categories.add(node['data']['label'])
    return node_categories

node_category_colors = {'Protein':'#a3b5c7',\
                        'BioProcess':'#4d2f9e',\
                            'MolFunc':'#7fe63a',\
                            'CellComp':'#f1c2d4',\
                            'Pathway':'#14ad9b',\
                            'Complex':'#ab3fd5',\
                            'Tissue':'#ef4b2c',\
                            'Domain':'#1e6723',\
                            'HomologousSuperFam':'#5638e7',\
                            'KinaseFam':'#94cd85',\
                            'KinaseFamGroup':'#af853e'}

def get_graph_style(elements):
    edge_labels = get_edge_labels(elements['edges'])
    node_labels = get_node_labels(elements['nodes'])
    node_styles = [NodeStyle(node_label, node_category_colors.get(node_label),'name') for node_label in node_labels]
    edge_styles = [EdgeStyle(edge_label, caption='label', directed=True) for edge_label in edge_labels]
    return node_styles, edge_styles

def build_inference_grid(data):
    data.columns = ['Kinase Gene','Kinase UniProt ID','Sites predicted (count)', 'p-value', 'Adjusted p-value']
    gb = GridOptionsBuilder.from_dataframe(data)
    gb.configure_pagination(enabled=True, paginationAutoPageSize=False, paginationPageSize=10)

    AgGrid(
            data, 
            gridOptions=gb.build(), 
            autoHeight=True,
            width='100%',
            allow_unsafe_jscode=True, #Set it to True to allow jsfunction to be injected
            enable_enterprise_modules=False, 
            fit_columns_on_grid_load=True,
            reload_data=False
            )