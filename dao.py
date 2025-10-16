import pandas as pd
import streamlit as st

from neo4j import GraphDatabase
import os

class DatabaseAccess:

    '''
    This class establishes connection to the database and provides a driver instance.
    '''
    neo4j_driver = None

    @classmethod
    def get_driver(cls):
        '''
        This method creates a new database driver and verifies connectivity, if one doesn't exist. Else, it returns the existing driver object.
        '''
        if cls.neo4j_driver is None:
            cls.neo4j_driver = GraphDatabase.driver(os.getenv("DB_URI"), 
                                                    auth=(os.getenv("DB_USER"),os.getenv("DB_PWD")), 
                                                    database=os.getenv("DB_NAME"))
            cls.neo4j_driver.verify_connectivity()
        return cls.neo4j_driver
    
    @classmethod
    def close_driver(cls):
        '''
        This method closes the existing driver object
        '''
        if cls.neo4j_driver:
            cls.neo4j_driver.close()
            cls.neo4j_driver = None


class KinaseSubstrateDAO:

    QUERY_MAX_IDT_DICT = {1:"""
                            MATCH path=(k:Protein)-[k_r]->(idt)<-[s_r]-(s:Protein) 
                            WHERE k.id='{k_id}' and s.id='{s_id}' 
                            RETURN k.name as k_name, k.id as k_id, labels(k) as k_label,
                                type(k_r) as k_r, elementId(k_r) as k_r_id,
                                idt.name as idt_name, idt.id as idt_id, labels(idt) as idt_label,
                                type(s_r) as s_r, elementId(s_r) as s_r_id,
                                s.name as s_name, s.id as s_id,labels(s) as s_label,
                                '' as idt1_r, '' as idt1_r_id,
                                '' as idt1_name, '' as idt1_id,"" as idt1_label,
                                '' as idt_type
                            """,
                          2:"""
                            MATCH(k:Protein)-[k_r1]->(k_idt)-[k_r]->(idt)<-[s_r]-(s:Protein) 
                            WHERE k.id='{k_id}' and s.id='{s_id}' 
                            RETURN k.name as k_name, k.id as k_id,labels(k) as k_label,
                                type(k_r) as k_r, id(k_r) as k_r_id,
                                idt.name as idt_name, idt.id as idt_id,labels(idt) as idt_label,
                                type(s_r) as s_r, id(s_r) as s_r_id,
                                s.name as s_name, s.id as s_id,labels(s) as s_label,                        
                                type(k_r1) as idt1_r, id(k_r1) as idt1_r_id,
                                k_idt.name as idt1_name, k_idt.id as idt1_id,labels(idt1) as idt1_label,
                                'k_idt1' as idt_type
                            UNION
                            MATCH(k:Protein)-[k_r]->(idt)<-[s_r]-(s_idt)<-[s_r1]-(s:Protein) 
                            WHERE k.id='{k_id}' and s.id='{s_id}'
                            RETURN k.name as k_name, k.id as k_id,labels(k) as k_label,
                                type(k_r) as k_r,id(k_r) as k_r_id,
                                idt.name as idt_name,idt.id as idt_id,labels(idt) as idt_label,
                                type(s_r) as s_r,id(s_r) as s_r_id,
                                s.name as s_name,s.id as s_id,labels(s) as s_label,    
                                type(s_r1) as idt1_r,id(s_r1) as idt1_r_id,
                                s_idt.name as idt1_name,s_idt.id as idt1_id,labels(idt1) as idt1_label,
                                's_idt1' as idt_type
                            """}

    def __init__(self, db_driver) -> None:
        '''
        This is the initialization method for KinaseSubstrateDAO. It requires a database driver object to be instantiated.
        
        Attributes:
        db_driver: neo4j driver object
        '''
        self.driver = db_driver

    def query_triples(self,tx,query):
        triples = set()
        for record in tx.run(query):
            record_parser = DBRecordParser(record)
            triples.update(record_parser.get_triples())
        return triples     

    def get_ks_links(self, kinase,substrate,intermediate_max=1):
        query = KinaseSubstrateDAO.QUERY_MAX_IDT_DICT[intermediate_max]
        query = query.format(k_id=kinase,s_id=substrate)
        ks_triples = set()
        with self.driver.session() as session:
            ks_triples = session.execute_read(self.query_triples, query)
        return ks_triples

class DBRecordParser:

    def __init__(self,record):
        self.k_name,self.k_id,self.k_category = record['k_name'],record['k_id'],record['k_label'][0]
        self.k_r,self.k_r_id = record['k_r'],record['k_r_id']
        self.idt_name,self.idt_id,self.idt_category = record['idt_name'],record['idt_id'],record['idt_label'][0]
        self.s_r,self.s_r_id = record['s_r'],record['s_r_id']
        self.s_name,self.s_id,self.s_category = record['s_name'],record['s_id'],record['s_label'][0]
        self.idt1_r,self.idt1_r_id = record['idt1_r'],record['idt1_r_id']
        self.idt1_name,self.idt1_id,self.idt1_category = record['idt1_name'],record['idt1_id'],record['idt1_label']
        self.idt_type = record['idt_type']
        self.triples = set()

    def _parse_to_triples(self):

        triples = set()
        k_node = Node(self.k_id,self.k_name,self.k_category)
        idt_node = Node(self.idt_id,self.idt_name,self.idt_category)
        s_node = Node(self.s_id,self.s_name,self.s_category)

        k_rel = Edge(self.k_r_id,self.k_r)
        s_rel = Edge(self.s_r_id,self.s_r)
        if self.idt_type:
            if self.idt_type == 'k_idt1':
                k_idt_node = Node(self.idt1_id,self.idt1_name,self.idt1_category[0])
                idt_rel = Edge(self.idt1_r_id,self.idt1_r)
                k_to_k_idt = Triple(k_node,k_rel,k_idt_node)
                k_idt_to_idt = Triple(k_idt_node,idt_rel,idt_node)
                s_to_idt = Triple(s_node,s_rel,idt_node)
                triples.add(k_to_k_idt)
                triples.add(k_idt_to_idt)
                triples.add(s_to_idt)
            elif self.idt_type == 's_idt1':
                s_idt_node = Node(self.idt1_id,self.idt1_name,self.idt1_category[0])
                idt_rel = Edge(self.idt1_r_id,self.idt1_r)
                k_to_idt = Triple(k_node,k_rel,idt_node)
                s_to_s_idt = Triple(s_node,s_rel,s_idt_node)
                s_idt_to_idt = Triple(s_idt_node,s_rel,idt_node)
                triples.add(k_to_idt)
                triples.add(s_to_s_idt)
                triples.add(s_idt_to_idt)
        else:
            k_to_idt = Triple(k_node,k_rel,idt_node)
            s_to_idt = Triple(s_node,s_rel,idt_node)
            triples.add(k_to_idt)
            triples.add(s_to_idt)
        return triples
    
    def get_triples(self):
        if not self.triples:
            self.triples = self._parse_to_triples()
        return self.triples


class Triple:

    def __init__(self, head, edge, tail):
        self.head = head
        self.rel = edge
        self.tail = tail

    def get(self):
        return self.head, self.rel, self.tail
    
    def __eq__(self,other):
        if isinstance(other,Triple):
            return (self.head == other.head)  \
                    & (self.tail == other.tail) \
                    & (self.rel == other.rel)
        return False

    def __hash__(self):
        return hash((self.head,self.rel,self.tail))
    
class Node:

    def __init__(self,id,name,category):
        self.id = id
        self.name = name
        self.category = category

    def get_id(self):
        return self.id
    
    def get_name(self):
        return self.name
    
    def get_category(self):
        return self.category
    
    def __eq__(self,other):
        if isinstance(other,Node):
            return self.id == other.id
        return False
    
    def __hash__(self):
        return hash(self.id)

    
class Edge:

    def __init__(self,id,label):
        self.id = id
        self.label = label

    def get_id(self):
        return self.id
    
    def get_label(self):
        return self.label
    
    def __eq__(self,other):
        if isinstance(other,Edge):
            return self.id == other.id
        return False
    
    def __hash__(self):
        return hash(self.id)


# '''
# This method returns all kg substrates as a list of tuples. Each tuple represents a substrate protein, the phosphorylation site and the -/+4 motif
# '''
# @st.cache_data
# def get_kg_substrates():
#     kg_substrates_df = pd.read_csv('data/substrates_motif.csv', sep='|')
#     substrate_proteins,sites,motifs = kg_substrates_df['Seq_Substrate'].to_list(), kg_substrates_df['Site'].to_list(), kg_substrates_df['Motif'].to_list()
#     return [(substrate_proteins[idx],sites[idx],motifs[idx]) for idx in range(0,kg_substrates_df.shape[0])]

'''
This method returns all kg substrates as a list of tuples. Each tuple represents a substrate protein, the phosphorylation site and the -/+4 motif
'''
@st.cache_data
def get_kg_substrates():
    kg_substrates_df = pd.read_csv('data/substrates_motif.csv', sep='|')
    substrate_genes, substrate_proteins,sites,motifs = kg_substrates_df['gene'].to_list(), kg_substrates_df['Seq_Substrate'].to_list(), kg_substrates_df['Site'].to_list(), kg_substrates_df['Motif'].to_list()
    return [(substrate_genes[idx],substrate_proteins[idx],sites[idx],motifs[idx]) for idx in range(0,kg_substrates_df.shape[0])]
