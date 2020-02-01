from os.path import commonprefix
from FAdo.reex import *
from FAdo.cfg import *
from FAdo.fa import *

# Tree Structure
class eTree:
    def __init__(self,regex):
        ReExp = concat(atom("#"),concat(regex.marked(),atom("$")))
        self.NT = dict()
        self.atoms = set()
        self.next = dict()
        self.sub_trees = dict()
        self.alpha = set()
        self.full_tree = construct(ReExp,"1",self)
        self.root = self.NT.get("121")
        self.Color()
    
    def followList(self):
        l = dict()
        for x in self.atoms:
            f_list = set()
            for w in self.atoms:
                if w.type!=epsilon:
                    if w.follow(x):
                        f_list.add(w)
            l[x]=f_list
        return l

    def build_ta(self,a):
        if self.sub_trees.get(a) is not None:
            return self.sub_trees.get(a)
        ta_atoms = set()
        for n in self.atoms:
            if n.exp.val[0]==a:
                ta_atoms.add(n)
        ta = set()
        for x in ta_atoms:
            for y in ta_atoms:
                if x is y:
                    ta.add(x.parent)
                else:
                    ta.add(self.NT.get(commonprefix([x.path,y.path])))
        self.sub_trees[a] = ta_atoms.union(ta.union(self.find_stars(a)))
        return self.sub_trees.get(a)

    def find_stars(self,a):
        ta = set()
        for n in self.atoms:
            if n.exp.val[0] == a:
                n = n.parent
                while n != self.root:
                    if n.type is star:
                        ta.add(n)
                    n = n.parent
        return ta
    
    def has_Next(self,n,a):
        s = self.next.get(n)
        if s is None:
            return False
        for x in s:
            if x.exp.val[0] == a:
                return True
        return False
    
    def Color(self):
        for x in self.atoms:
            x.pSupFirst().parent.colors.add(x.exp.val[0])
        return
    
    def isDeterministic(self):
        for x in self.alpha:
            buildNext(x,self.root,set())
        if not self.cond_P1() or not self.cond_P2():
            return False
        return self.check_determinism(self.root)
    
    def cond_P1(self):
        for a in self.atoms:
            for b in self.atoms:
                if a is b:
                    continue
                if a.pSupFirst() == b.pSupFirst():
                    if a.exp.val[1] == b.exp.val[1]:
                        print "P1"
                        return False
        return True
    
    def cond_P2(self):
        for a in self.alpha:
            ta = self.build_ta(a)
            for n in ta:
                k = n.FollowAfter
                count = 0
                for x in k:
                    if x.exp.val[0] == a:
                        count = count + 1
                if count > 1:
                    print "P2"
                    return False
        return True
    
    def check_determinism(self,n):
        for x in n.colors:
            if not checkNode(n,x):
                return False
        if n.left is not None:
            l = self.check_determinism(n.left)
        else:
            return True
        if n.right is not None:
            return (l and self.check_determinism(n.right))
        return l

# Searches for non-determinism in each node
def checkNode(n,a):
    te = n.tree
    ta = te.build_ta(a)
    f = n.FirstPos(ta)
    s = n.LSA()
    if n.right.exp.ewp():
        if te.has_Next(n,a):
            return False
        if s is None:
            return True
        if s.FirstPos(ta) is f and s.reflexive(n.pSupLas()):
            return False
    return True
                 
# Constructor for eTree structure. Construct nodes.
def construct(reg_exp,path,t):
        node = Node(reg_exp,path,t)
        t.NT[node.path] = node
        if node.type is star:
            node.left = construct(node.exp.arg,path+'1',t)
        elif node.type is concat or node.type is disj:
            node.left = construct(node.exp.arg1,path+'1',t)
            node.right = construct(node.exp.arg2,path+'2',t)
        elif node.type is position:
            sym = node.exp.val[0]
            t.alpha.add(sym)
        return node

class Node:
    def __init__(self,expr,string,t):
        self.type  = type(expr)
        self.path  = string
        self.exp   = expr
        self.tree  = t
        self.right = None
        self.left  = None
        self.first = set()
        self.last  = set()
        self.lsa   = None
        self.supfirst = None
        self.suplast = None
        self.FollowAfter = set()
        self.colors = set()
        if self.type is position:
            t.atoms.add(self)
        if string == '1':
            self.parent = None
        else:
            self.parent = t.NT.get(string[:-1])
    
    # check if n is ancestor of self
    def reflexive(self,n):
        if n is None:
            return False
        if commonprefix([self.path,n.path])==n.path: #LCA(self,n)
            return True
        return False
        
    # boolean that states if an atom is the first from the right child of a concatenation
    def SupFirst(self):
        father = self.parent
        if father is not None and father.type is concat:
            if father.right is self and not father.left.exp.ewp():
                return True
        return False

    # boolean that states if an atom is the last from the left child of a concatenation
    def SupLast(self):
        father = self.parent
        if father is not None and father.type is concat:
            if father.left is self and not father.right.exp.ewp():
                return True
        return False
    
    # pointer to the first(right child) of a concatenation
    def pSupFirst(self):
        if self.supfirst is not None:
            return self.supfirst
        self.supfirst = None
        if self.SupFirst() or self.path=="1":
            self.supfirst = self
        else:
            self.supfirst = self.parent.pSupFirst()
        return self.supfirst

    # pointer to the last(left child) of a concatenation
    def pSupLast(self):
        if self.suplast is not None:
            return self.suplast
        if self.SupLast() or self.path=="1":
            self.suplast = self
        else:
            self.suplast = self.parent.pSupLast()
        return self.suplast
    
    #commonprefix is the Python prelude function used for LCA
    def follow(self,p):
        node1 = self.tree.NT.get(commonprefix([p.path,self.path]))
        if node1 is None:
            return False
        if node1.type is concat:
            if self.follow_concat(p,node1):
                return True
        node = node1.LSA()
        if node is not None:
            if self.follow_star(p,node):
                return True
        return False
    
    # Follow by Concatenation
    def follow_concat(self,p,lca):
        if self in lca.right.First() and p in lca.left.Last():
            return True
        return False
    
    # Follow by Star           
    def follow_star(self,p,lsa):
        if self in lsa.First() and p in lsa.Last():
            return True
        return False
    
    def First(self):
        if self.first != set():
            return self.first
        for node in self.tree.atoms:
            if isFirst(self,node):
                self.first.add(node)
        return self.first
    
    def Last(self):
        if self.last != set():
            return self.last
        for node in self.tree.atoms:
            if isLast(self,node):
                self.last.add(node)
        return self.last
    
    def LSA(self):
        if self.lsa is not None:
            return self.lsa
        if self.type is star:
            self.lsa = self
        elif len(self.path) == 1:
            self.lsa = None
        else:
            self.lsa = self.parent.LSA()
        return self.lsa
    
    def FirstPos(self,ta):
        k = set()
        for x in ta:
            if x in self.First():
                k.add(x)
        return k
    
    def follow_after(self):
        s = set()
        for p in self.Last():
            for q in self.tree.atoms:
                if q.follow(p) and not q.reflexive(self):
                    s.add(q)
        self.FollowAfter = s
        return s

def isFirst(n,p):
    k = p.pSupFirst()
    if k is not None:
        if p.reflexive(n) and n.reflexive(k):
            return True
    return False

def isLast(n,p):
    k = p.pSupLast()
    if k is not None:
        if p.reflexive(n) and n.reflexive(k):
            return True
    return False

def buildNext(a,n,y):
        if n.SupLast():
            y = set()
        te = n.tree
        ta = n.tree.build_ta(a)
        if ta is None:
            return False
        dad = get_parent(n,ta)
        if dad is not None:
            r = get_right(dad,ta)
            if dad.type is concat and get_left(dad,ta) is n and r is not None:
                if not n.SupLast() or get_parent(n,ta) is n.parent:
                    y = y.union(r.FirstPos(ta))
        n.tree.next[n] = set()
        for x in y:
            if not x.reflexive(n):
                n.tree.next[n].add(x)
        if n.type is star:
            y = y.union(n.FirstPos(ta))
        if len(y)>2:
            return False
        l = get_left(n,ta)
        if l is None or l.exp.ewp():
            return True
        else:
            b = buildNext(a,l,y)
        r = get_right(n,ta)
        if r is None or r.exp.ewp():
            return b
        return (b and (buildNext(a,r,y)))

def get_parent(n,tree):
    t = n.tree
    while n is not t.root:
        parent = n.parent
        if parent in tree:
            return parent
        n = n.parent
    return None

def get_left(n,tree):
    if n.type is position:
        return None
    path = n.path+'1'
    possible = set()
    for x in tree:
        if len(x.path)>=len(path) and commonprefix([x.path,path])==path:
            possible.add(x.path)
    if len(possible) == 0:
        return None
    return n.tree.NT.get(min(possible))

def get_right(n,tree):
    if n.type is position:
        return None
    path = n.path+'2'
    possible = set()
    for x in tree:
        if len(x.path)>=len(path) and commonprefix([x.path,path])==path:
            possible.add(x.path)
    if len(possible) == 0:
        return None
    return n.tree.NT.get(min(possible))