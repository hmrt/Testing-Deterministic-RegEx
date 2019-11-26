from FAdo.reex import *
import re

NT = dict({})
class Node:
    
    def __init__(self,expr,string):
        self.type = type(expr)
        self.path = string
        self.exp = expr
        NT[self.path] = self
    
    def set_left(self,n):
        NT[self.path+'1'] = n
        
    def set_right(self,n):
        NT[self.path+'2'] = n

    def get_left(self):
        return NT.get(self.path+'1')
    
    def get_right(self):
        return NT.get(self.path+'2')
        
    def get_parent(self):
        return NT[self.path[:-1]]
    
    def printTree(self):
        if NT.get(self.path):
            print self.exp
            if NT.get(self.path+'1'):
                l = NT[self.path+'1']
                l.printTree()
            if NT.get(self.path+'2'):
                r = NT[self.path+'2']
                r.printTree()
        
    def get_root(self):
        return NT.get("1")
    
    def check_reex(self):
        if self.type == 'atom':
            if self.get_left() or self.get_right():
                return False
        elif self.type == 'star':
            if self.get_left() is None or self.get_right():
                return False
            if self.get_left().check_reex() is False:
                return False
        else:
            if self.get_left() is None or self.get_right() is None:
                return False
            if self.get_left().check_reex() is False or self.get_right().check_reex() is False:
                return False
        return True

def get_root():
    return NT.get("1")

def LCA(p,q):
    str1 = p.path
    str2 = q.path
    ret = ""
    i = 0
    while i < len(str1) and i < len(str2):
        if str1[i]!=str2[i]:
            break
        ret = ret+str1[i]
        i = i+1
    return ret

def eTree(reg_exp,path):
    node = Node(reg_exp,path)
    if type(reg_exp)==atom:
        return
    elif type(reg_exp)==star:
        son = read_reex(node.exp.arg,path+"1")
    else:
        son1 = read_reex(node.exp.arg1,path+"1")
        son2 = read_reex(node.exp.arg2,path+"2")
				
def lemma2(p,q):
    node1 = NT.get(LCA(p,q))
    if type(node1)==concat and q in node1.get_right().first() and p in node1.get_left().last():
        return True
    while type(node1)!=star:
        if node1.path=="1":
            return False
        node1 = node1.get_parent()
    if q in node1.first() and p in node1.last():
        return True
    return False
