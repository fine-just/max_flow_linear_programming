import string
from numpy import random
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import itertools
#παιρνει τα κλειδια απο λεξιλόγιο
def get_keys(dict1):
    keys=[]
    
    for i in dict1:
        if(dict1[i]!=0):
            keys.append(i)
        
    return keys
#βρίσκει όλες τις διαδρομές απο δίκτυο
#που έχουν θετικές χωρητικότητες
def all_paths(graph,start,end,path,paths):
    #μολις φτάσει στην τελική κορυφή
    #επιστρέφει την τελική κορυφή
    if(start==end):
        return end
    else:
        #ανανεώνει το path
        path=path+[start]
        #ελέγχει αν έχει παιδιά η κορυφή
        check_key=start in graph.keys()
        if(check_key==False):
            return None
            
        child=get_keys(graph[start])
        num=0
        for i in child:
            #αντιγράφει την διαδρομή
            #για να επαναληφθεί για τα
            #παιδιά του start
            path1=path.copy()
            #αμα πάει να γίνει loop κορυφών
            #επιστρέφει None
            if(i in path1):
                return None
            else:
                #κάνει αναδρομή με το παιδί του start
                test=all_paths(graph,i,end,path1,paths)
                if(test!=None):
                
                    path1.append(test)
                    paths.append(path1)
            #print(paths)
            #μόλις βρεθούν 2 διαδρομές
            #τελειώνει η διαδικασία
            if(len(paths)==2):
                break
                
            
#------start of class graphs----------------
class graphs:
    #ορίζουμε τα αρχικά χαρακτηριστικά
    def __init__(self):
        self.forward={}
        self.backward={}
        self.max_flow_ff=0
        self.source=''
        self.sink=''
        self.paths=[]
        self.variables=[]
        self.equations={}
        self.anisosi={}
    #συνάρτηση που κάνει την μέθοδο Simplex
    #για ένα πίνακα r5
    def simplex(self,r5):
        
        r6=r5.copy()
        #βρίσκει τα όρια του πίνακα
        fin_x=len(r5[0,:])-1
        fin_y=len(r5[:,0])-1

        #briskw arn sto r5[fin_y,:fin_x]
        #μέχρι να γίνει break
        while(True):
            #df=pd.DataFrame({'x1':r5[:,0],'x2':r5[:,1],'x3':r5[:,2],'x4':r5[:,3],'x5':r5[:,4],'x6':r5[:,5],'x7':r5[:,6],'y1':r5[:,7],'y2':r5[:,8],'y3':r5[:,9],'y4':r5[:,10],'y5':r5[:,11],'y6':r5[:,12],'Z':r5[:,13],'b':r5[:,14]})
            #df1=pd.DataFrame({'x1':r5[:,0],'x2':r5[:,1],'x3':r5[:,2],'x4':r5[:,3],'x5':r5[:,4],'x6':r5[:,5],'x7':r5[:,6],'b':r5[:,14]})

            #βρίσκει το ελάχιστο απ'τη
            #τελευταία στήλη
            ind1=np.argmin(r5[fin_y,:fin_x])
            #print(ind1)
            #αν είναι θετικό
            #τελειώνει η διαδικασία
            if(r5[fin_y,ind1]>=0):
                break
            #βρίσκει την στήλη του ind1
            min_stil=r5[:fin_y,ind1]
            base=r5[:fin_y,fin_x]
            b2=base.copy()
            min1=1000
            min_ind=-1
            #βρίσκει το ελάχιστο της στήλης
            for i in range(len(min_stil)):
                if(min_stil[i]!=0):
                    g=b2[i]/min_stil[i]
                    if(min1>g and g>0):
                        min1=g
                        min_ind=i
            #print(min1)          
            #ato min_ind tha ginei pivot           
            r5[min_ind]=r5[min_ind]/r5[min_ind,ind1]
            for i in range(fin_y+1):
                if(i!=min_ind):
                    r5[i]=r5[i]-r5[min_ind]*r5[i,ind1]
            #print(df)
            #print('--------------')
        #r5=np.round(r5, 3)
        #df=pd.DataFrame({'x1':r5[:,0],'x2':r5[:,1],'x3':r5[:,2],'x4':r5[:,3],'y1':r5[:,4],'y2':r5[:,5],'y3':r5[:,6],'y4':r5[:,7],'Z':r5[:,8],'b':r5[:,9]})
        #print('g1='+str(g1)+'g2='+str(g2))
        #print(df)
        
        return r5
    #ορίζουμε τις ακμές του δικτύου
    def insert_edges(self,start,end,cost):
        #εισάγει τις κορυφές και την χωρητικότητα
        #στο self.forward
        check_f=start in self.forward.keys()
        check_b=end in self.backward.keys()
        if(check_f==False):
            self.forward[start]={}
        self.forward[start][end]=cost

        if(check_b==False):
            self.backward[end]={}
        self.backward[end][start]=0
    
    def get_dict_keys(self,dict1):
        keys=[]
    
        for i in dict1:
            
            keys.append(i)
        
        return keys
    def forward_cost(self,start,end,cost):
        check_f=start in self.forward.keys()
        if(check_f==False):
            return 'the start edge doesnt exists'
        else:
            check_f1=end in self.forward[start].keys()
            if(check_f1==False):
                return 'the end edge doesnt exists'
            self.forward[start][end]=cost
            
        
    def backward_cost(self,start,end,cost):
        check_f=start in self.backward.keys()
        if(check_f==False):
            return 'the start edge doesnt exists'
        else:
            check_f1=end in self.backward[start].keys()
            if(check_f1==False):
                return 'the end edge doesnt exists'
            self.backward[start][end]=cost
    def insert_source_sink(self,start,end):
        self.source=start
        self.sink=end

    def find_paths(self):
        self.paths=[]
        path=[]
        paths=self.paths
        a=self.forward
        start1=self.source
        end1=self.sink
        all_paths(a,start1,end1,path,paths)


    def find_forward_cost(self,start,end):
        check_f=start in self.forward.keys()
        if(check_f==False):
            return None
        else:
            check_f1=end in self.forward[start].keys()
            if(check_f1==False):
                return None

        return p1.forward[start][end]
    def find_backward_cost(self,start,end):
        check_f=start in self.backward.keys()
        if(check_f==False):
            return None
        else:
            check_f1=end in self.backward[start].keys()
            if(check_f1==False):
                return None
        return p1.backward[start][end]
    def insert_graph(self):
        stop=False
        print('input data as the example')
        print('nodes from to :T,S')
        print('with cost:45')
        while(stop==False):
            
            print('nodes from to:')
            nodes=list(input())
            nodes.remove(',')
            print('/n')
            print('with cost:')
            cost=int(input())
            self.insert_edges(nodes[0],nodes[1],cost)
            print('/n')
            print('Do you want to stop? Y/N \n')
            stop1=input()
            if(stop1=='Y'):
                stop=True
        print('the graph starts from point:')
        start=input()
        
        print('/n')
        print('the graph ends in point:')
        end=input()
        
        self.insert_source_sink(start,end)

    def instert_array_graph(self,array,source,sink):
        for i in array:
            self.insert_edges(i[0],i[1],i[2])
        self.insert_source_sink(source,sink)
        
#ford and fulkerson method----------------------------------------------------
    def f_and_f_method(self):
        start = time.time()
        #βρίσκω τα 2 πρώτα paths
        self.find_paths()
        
        iterations=1
        #μέχρι να μην υπάρχουν άλλα paths
        #στο δίκτυο 
        while(self.paths!=[]):
            augmentations=0
            #παίρνει το πρώτο path
            curr_path=self.paths[0]
            #αρχικοποιεί τον πίνακα curr_cost
            curr_cost=[]
            #αρχικοποιεί το min_cost
            min_cost=0
            #αποθηκεύει όλες τις χωρητικότητες
            #της διαδρομής στο curr_cost
            for i in range(len(curr_path)-1):
                curr_cost.append(self.forward[curr_path[i]][curr_path[i+1]])
            #βρίσκει την ελάχιστη 
            min_cost=min(curr_cost)
            #ανανεώνει τα forward flow και τα backward flow
            augmentations=augmentations+1
            for j in range(len(curr_path)-1):
                
                self.forward[curr_path[j]][curr_path[j+1]]=self.forward[curr_path[j]][curr_path[j+1]]-min_cost
                
                self.backward[curr_path[j+1]][curr_path[j]]=self.backward[curr_path[j+1]][curr_path[j]]+min_cost
           
            print('at iteration %d we have %d augmentations'%(iterations,augmentations))
            self.find_paths()
            iterations=iterations+1
        m=0
        #βρίσκει το max_flow
        for k in self.backward[self.sink]:
            m=m+self.backward[self.sink][k]

        self.max_flow_ff=m
        print('augmentations',augmentations)
        print('iterations',iterations)
        end = time.time()
        print(f"Runtime of the Ford and Fulkerson method is {end - start}")
        return m
#FOR SIMPLEX METHOD----------------------------------------------------------
    def find_all_variables(self):
        for i in self.forward:
            if((i in self.variables)==False):
                self.variables.append(i)
        for j in self.backward:
            if((j in self.variables)==False):
               
                self.variables.append(j)
    def find_equations(self):
        self.find_all_variables()
        for i in self.variables:
            self.equations[i]={'1':[],'-1':[]}

            if((i in self.backward)==True):
                
                #χρησιμοποιοιυμε την self.get_dict_keys
                #γτ η def get_keys παιρνει τα κλειδια μονο για κοστος!=0 
                keys_b=self.get_dict_keys(self.backward[i])
                
                for k in keys_b:
                    self.equations[i]['-1'].append(k)

            if((i in self.forward)==True):
                keys_f=self.get_dict_keys(self.forward[i])
                
                for j in keys_f:
                    self.equations[i]['1'].append(j)
        self.equations[self.source]['-1'].append(self.sink)
        self.equations[self.sink]['1'].append(self.source)

    def simplex_method(self):
        start=time.time()
        self.find_equations()
        self.anisosi['m']={}
        for i in self.equations:
            
            if(i!=self.sink and i !=self.source):
                #print(i)
                positive_keys=self.equations[i]['1']
                negative_keys=self.equations[i]['-1']
                pos_min=0
                neg_min=0
                more_keys_p=0
                pos_keys=[]
                neg_keys=[]
                more_keys_n=0
                
                if(len(positive_keys)>1):
                    more_keys_p=1
                if(len(negative_keys)>1):
                    more_keys_n=1


                for j in positive_keys:
                    pos_keys.append(i+j)
                        #print(i+j)
                    
                    pos_min=pos_min+self.forward[i][j]
                j=0
                for j in negative_keys:
                    neg_keys.append(j+i)
                    
                    neg_min=neg_min+self.forward[j][i]
                
                if(neg_min<pos_min):
                    min1=neg_min
                else:
                    min1=pos_min

                for j in positive_keys:
                    if(self.forward[i][j]<min1):
                        self.anisosi[i+j]=self.forward[i][j]
                    else:
                        self.anisosi[i+j]=min1
                for j in negative_keys:
                    if(self.forward[j][i]<min1):
                        self.anisosi[j+i]=self.forward[j][i]
                    else:
                         self.anisosi[j+i]=min1
                if(more_keys_p==1):
                    str1=''
                    self.anisosi['m'][str(pos_keys)]=pos_keys,min1
                if(more_keys_n==1):
                    str2=''
                    self.anisosi['m'][str(neg_keys)]=neg_keys,min1

        #παιρνω τα κλειδια απο το source
        sink_source=[]
        keys_source=self.get_dict_keys(self.forward[self.source])
        source_min=0
        sink_min=0
        sour=[]
        sink=[]
        for i in keys_source:
            sink_source.append(self.source+i)
            sour.append(self.source+i)
            source_min=source_min+self.forward[self.source][i]
        keys_sink=self.get_dict_keys(self.backward[self.sink])
        for i in keys_sink:
            sink_min=sink_min+self.forward[i][self.sink]
            sink.append(i+self.sink)
            #sink_source.append(i+self.sink)
        min_s=min(sink_min,source_min)
        
        self.anisosi['m'][str(sink)]=sink,min_s
        self.anisosi['m'][str(sour)]=sour,min_s

        
        #στο self.anisosi ειναι ολες οι ανισωσεις
        #θα φτιαξω τους πινακες
        #βρισκω τον αριθμο των μεταβλητων
        agnwstoi=len(self.anisosi)-1
        #συνολο ανισωσεων
        num_aniso=agnwstoi+len(self.anisosi['m'])
        I_pinakas=np.eye(num_aniso+1)
        r_array=[]

        #θα εχω ολους τους κομβους σε ενα πινακα
        
        node_keys=self.get_dict_keys(self.anisosi)
        node_keys.remove('m')
        nodes={}
        for i in range(len(node_keys)):
            nodes[node_keys[i]]=i
        

        array_aniso=np.zeros((num_aniso+1,agnwstoi))
        apotel=np.zeros((num_aniso+1,1))
        row=0
        for i in self.anisosi:
            if(i=='m'):
                for j in self.anisosi[i]:
                    for k in self.anisosi[i][j][0]:
                        array_aniso[row,nodes[k]]=1
                    apotel[row]=self.anisosi[i][j][1]
                    row=row+1
            else:
                array_aniso[row,nodes[i]]=1
                apotel[row]=self.anisosi[i]
                row=row+1
        
        for i in sink_source:
            #array_aniso[-1,nodes[i]]=-0.5
            array_aniso[-1,nodes[i]]=-1.0

        #φτιαχνω τον τελικο πινακα
        f1=np.concatenate((array_aniso,I_pinakas),axis=1)
        final_ar=np.append(f1,apotel,axis=1)
        #return final_ar
        #θα τρεξει με την simplex
        self.simplex(final_ar)
        data={}
        for i in nodes:
            data[i]=final_ar[:,nodes[i]]
        data['TS']=final_ar[:,len(final_ar[0,:])-2]
        data['b']=final_ar[:,-1]
        self.df_simplex=pd.DataFrame(data)
        self.max_flow_simplex=(self.df_simplex['b'][len(self.df_simplex['b'])-1])
        end = time.time()
        print(f"Runtime of the simplex method  is {end - start}")
        return self.max_flow_simplex
#for Edmond karp method-----------------------------------------------------------------------
    def bfs_search(self,graph):
        start=self.source
        end=self.sink
        visited=[]
        queue=[[start]]
        while queue:
            path=queue.pop(0)
            node=path[-1]
            if((node in visited)==False):
                child=graph[node]
                
                for j in child:
                    new_path=list(path)
                    new_path.append(j)
                    queue.append(new_path)
                    if(j==end):
                        return new_path
                visited.append(node)

    def bfs_second_version(self):
        start1=time.time()
        start=self.source
        end=self.sink
        visited=[]
        queue=[[start]]
        while queue:
            path=queue.pop(0)
            
            node=path[-1]
            
            if((node in visited)==False and (node in self.forward)==True):
                
                child=get_keys(self.forward[node])
                
                for j in child:
                    new_path=list(path)
                    new_path.append(j)
                    queue.append(new_path)
                    if(j==end):
                        end1=time.time()
                        return new_path
                visited.append(node)
    
        
    def edmond_graph(self):
        graph={}
        
        self.find_all_variables()
        for i in self.variables:
            if((i in self.forward)==True):
                graph[i]=get_keys(self.forward[i])
            else:
                graph[i]=[]
        graph[self.sink]=[]
        return graph
    def edmond_karp_method(self):
        start=time.time()
        sh_path=self.bfs_second_version()
        iterations=1
        augmentations=0
        while(sh_path!=None):
        
            min1=[]
            for i in range(len(sh_path)-1):
                min1.append(self.forward[sh_path[i]][sh_path[i+1]])
            min_capacity=min(min1)
            augmentations=augmentations+1
            for i in range(len(sh_path)-1):
                self.forward[sh_path[i]][sh_path[i+1]]=self.forward[sh_path[i]][sh_path[i+1]]-min_capacity
                self.backward[sh_path[i+1]][sh_path[i]]=self.backward[sh_path[i+1]][sh_path[i]]+min_capacity
            print('at iteration %d we have %d augmentations'%(iterations,augmentations))
            sh_path=self.bfs_second_version()
            iterations=iterations+1
            augmentations=0
        max_flow=0
        for i in self.backward[self.sink]:
            max_flow=max_flow+self.backward[self.sink][i]
        end = time.time()
        print(f"Runtime of the program is {end - start}")
        print(max_flow)
        return max_flow
    
#AL-AMIN-KHAN METHOD-------------------------------------------------
    def all_the_paths(self,graph,start,end,path,paths):
        if(start==end):
            return end
        else:
            path=path+[start]
            check_key=start in graph.keys()
            if(check_key==False):
                return None
            
            
            child=get_keys(graph[start])
            num=0
            for i in child:
            
                path1=path.copy()
                if(i in path1 and path1.count(i)>=2):
                    
                    return None
                else:
                    
                    test=self.all_the_paths(graph,i,end,path1,paths)
                    if(test!=None):
                        
                
                        path1.append(test)
                        paths.append(path1)
    def paths_with_capacity(self):
        
        paths=[]
        path=[]
        
        a=self.forward
        start1=self.source
        end1=self.sink
        self.all_the_paths(self.forward,start1,end1,path,paths)
        path_capacity=[]
        final_path={}
        for curr_path in paths:
            path_capacity=[]
            for j in range(len(curr_path)-1):
                path_capacity.append(self.forward[curr_path[j]][curr_path[j+1]])
            min_cap=min(path_capacity)
            if min_cap in final_path:
                final_path[min_cap].append(curr_path)
            else:
                final_path[min_cap]=[curr_path]
        return final_path

    def find_al_path(self,paths,d):
        
        capacities=self.get_dict_keys(paths)
        capacities.sort(reverse=True)
        wanted_path=-10
        for i in capacities:
            if(i>=d):
                wanted_path=i
                break
        if(wanted_path==-10):
            return None

        else:
            #print(paths[wanted_path][0],capacities)
            return paths[wanted_path][0]

    def find_dict_sum(self,dict1):
        sum1=[]
        for i in dict1:
            sum1.append(dict1[i])
        return sum(sum1)
        

    def find_cap(self):
        comb=[]
        stuff=[]
        for i in self.forward:
            if(i not in stuff):
                stuff.append(i)
            k=p2.get_dict_keys(self.forward[i])
            for j in k:
                if(j not in stuff):
                    stuff.append(j)

        for L in range(0, len(stuff)+1):
            keys=[]
            for subset in itertools.combinations(stuff, L):
                keys=[]
                if (self.source in subset and self.sink not in subset):
                    for i in range(len(subset)-1):

                        if(subset[i] in p2.forward):
                            keys=keys+p2.get_dict_keys(self.forward[subset[i]])
                    if(subset[-1] in keys):
                        comb.append(subset)
        cap=[]
        for current_comb in comb:
            c1=0
            for i in current_comb:
                for j in self.forward[i]:
                    if(j not in current_comb):
                        c1=c1+self.forward[i][j]
            cap.append(c1)

        self.cap=cap.copy()
			    
			    
				    
        
            
    def al_amin_method(self):
        start1=time.time()
        self.find_all_variables()
        if(len(self.variables)<=20):
            self.find_cap()
        max_flow=0
        capacity=[]
        for i in self.forward:
            for j in self.forward[i]:
                capacity.append(self.forward[i][j])
        uc=max(capacity)
        lc=min(capacity)
        d=uc-lc
        end=0
        
        if(d==0):
            d=uc
        iterations=1
        augmentations=0
        while(end!=1):
            
            paths=self.paths_with_capacity()
            al_path=self.find_al_path(paths,d)
            
            if(al_path==None):
                while(True):
                    d=d/2
                    
                    paths=self.paths_with_capacity()
                    al_path=self.find_al_path(paths,d)
                    if(al_path!=None):
                        if(d<1):
                            
                            end=1
                        print('at iteration %d we have %d augmentations'%(iterations,augmentations))
                        iterations=iterations+1
                        augmentations=0
                        break
                    if(d<1):
                    
                        end=1
                        break


            current_cap=[]
            if(al_path!=None):
                for i in range(len(al_path)-1):
                    current_cap.append(self.forward[al_path[i]][al_path[i+1]])
                min_cap=min(current_cap)
                max_flow=max_flow+min_cap
                #print(max_flow)
                augmentations=augmentations+1
                for i in range(len(al_path)-1):
                    self.forward[al_path[i]][al_path[i+1]]=self.forward[al_path[i]][al_path[i+1]]-min_cap
                    self.backward[al_path[i+1]][al_path[i]]=self.backward[al_path[i+1]][al_path[i]]+min_cap
                
                
                print('at iteration %d we have %d augmentations'%(iterations,augmentations))
                iterations=iterations+1
                augmentations=0
            if(len(self.variables)<=20):
                if(max_flow in self.cap):
                    end=1
                    break
            
        print('finished')
        flow=[]
        for i in self.backward[self.sink]:
            flow.append(self.backward[self.sink][i])

        end1=time.time()
        print('time',end1-start1)
        print('max-flow',sum(flow))
        return max_flow
        #return sum(flow)
        

#-----------------end of class------------------------------------------

def make_graph(graph,numb_variables):
    #random.seed(1063974)
    upper=list(string.ascii_uppercase)
    lower=list(string.ascii_lowercase)
    upper.remove('T')
    upper.remove('S')
    numbers=list('0123456789')
    var=upper+lower+numbers
    final_var=[]
    nu=0
    if(numb_variables-2<=56):
        final_var=var.copy()
    else:
        while(True):
            v1=var[random.randint(len(var))]+var[random.randint(len(var))]
            if((v1 in final_var)==False):
                final_var.append(v1)
                nu=nu+1
            
            if(nu==numb_variables-2):
                break
  
    variables=final_var[:numb_variables-2]
    
    main=['T','S']
    connect_T=random.randint(1,numb_variables-2)
    connect_S=random.randint(1,numb_variables-2)
    i=0
    for i in range(connect_T):
        con=variables[random.randint(len(variables))]
        if((con in graph.backward)==False):
                
            graph.insert_edges('T',con,random.randint(1,20))
    
    for j in range(connect_S):
        con1=variables[random.randint(len(variables))]
            
        graph.insert_edges(con1,'S',random.randint(1,20))

    first_layer=get_keys(graph.forward['T'])
    new_var=variables.copy()
    for k in first_layer:
        new_var.remove(k)
    for k in first_layer:
        if(numb_variables==4):
            random_conn=1
        else:
            random_conn=random.randint(1,numb_variables-3)
        dont_touch={}
        if((k in graph.backward.keys())==True ):
            dont_touch=graph.backward[k]
        for k1 in range(random_conn):
            con2=new_var[random.randint(len(new_var))]
            if((con2 in dont_touch)==False ):
                
                graph.insert_edges(k,con2,random.randint(1,20))
                       
    graph.insert_source_sink('T','S')

def make_graph2(seed,graph,layers,numb_variables):
    
    numb=numb_variables
    l={}
    l1=round(numb_variables/layers)
    print(l1)
    
    connect_S=random.randint(1,numb_variables)

    variables=[]
    for  i in range(l1*layers):
        variables.append(str(i))
    #print(variables)

    
    
    

    
    new_var=variables.copy()

    laye={}
    b=0
    c=l1
    for i in range(layers):
        laye[i]=variables[b:c]
        b=b+l1
        c=c+l1
    for i in laye[0]:
        
        
         graph.insert_edges('T',i,random.randint(1,8))
         
    for j in range(connect_S):
        con1=variables[random.randint(len(variables))]
        if(con1 in graph.forward):
            if('S' not in graph.forward[con1]):
                graph.insert_edges(con1,'S',random.randint(1,20))
        else:
            graph.insert_edges(con1,'S',random.randint(1,20))
                
    for i in laye:
        for j in laye[i]:
            for k in range(random.randint(1,numb_variables-1)):
                
                con1=variables[random.randint(len(variables)-1)]
                if(con1!=j):
                    if(con1 in graph.backward):
                        if(j not in graph.backward):
                            graph.insert_edges(j,con1,random.randint(1,20))
                    else:
                        graph.insert_edges(j,con1,random.randint(1,20))
    graph.insert_source_sink('T','S')
                        
                    
                
                    
                    

                
    
    
            
        
def show_graph(graph):
    G = nx.DiGraph(directed=True)
  # Create empty graph
    graph.find_all_variables()
    G.add_nodes_from(graph.variables)  # Add nodes

    # Add edges\
    for i in graph.forward:
        for j in graph.forward[i]:
            G.add_edge(i,j, weight=graph.forward[i][j])
            
    
    
    
    # Create drawing
    
    pos = nx.spring_layout(G)  # List of positions of nodes
    weights = nx.get_edge_attributes(G, "weight") # List of weights
    nx.draw_networkx(G, pos, with_labels=True)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=weights)
    
    plt.title("Basic Graphs with Networkx")
    #plt.gcf().canvas.set_window_title("1")  # Hide window title

    # Display Graph
    plt.show()

def all_methods(first_arr,start,end):
    
    for i in range(5):
        first_level=graphs()
        #first_arr=[['S', '1', 18], ['S', '2', 5], ['3', '4', 4], ['S', '3', 2], ['1', '6', 20], ['1', '7', 10], ['2', '5', 20], ['4', '8', 30], ['5', '9', 5], ['5', '10', 4], ['6', '10', 12], ['7', '11', 11], ['7', '12', 13], ['8', '11', 14], ['9', '10', 6], ['10', '13', 7], ['10', '14', 14], ['11', '15', 9], ['11', '12', 8], ['12', '16', 15], ['13', 'T', 2], ['14', 'T', 19], ['15', '16', 11], ['16', 'T', 8]]
        first_level.instert_array_graph(first_arr,start,end)
        if(i==1):
            print('simplex')
            print(first_level.simplex_method())
            
            print('-----------\n')
            
        if(i==2):
            print('ford and fulkerson method')
            print(first_level.f_and_f_method())
            
            print('-----------\n')
        if(i==3):
            print('Edmond karp method')
            first_level.edmond_karp_method()
            
            print('-----------\n')
        if(i==4):
            print('Al-Amin Khan method')
            
            first_level.al_amin_method()
            
            print('-----------\n')
    first_level=graphs()
    first_level.instert_array_graph(first_arr,start,end)
    #show_graph(first_level)
   

if __name__ == "__main__":       
            
    p1=graphs()
    p1.insert_edges('T','A',10)


    p1.insert_edges('T','B',2)

    p1.insert_edges('A','C',5)
    p1.insert_edges('A','D',8)
    p1.insert_edges('B','D',7)
    p1.insert_edges('C','S',4)
    p1.insert_edges('C','A',4)
    p1.insert_edges('A','S',8)
    p1.insert_edges('D','S',1)
    p1.insert_edges('D','E',3)
    p1.insert_edges('E','S',11)
    p1.insert_source_sink('T','S')

    p2=graphs()
    p2.insert_edges('S','A',5)
    p2.insert_edges('S','B',6)
    p2.insert_edges('A','C',4)
    p2.insert_edges('A','D',2)
    p2.insert_edges('B','C',8)
    p2.insert_edges('D','T',4)
    p2.insert_edges('C','T',10)
            
    p2.insert_source_sink('S','T')



    g1=graphs()
    array_g1=[['0','1',8],['0','2',6],['1','2',2],['1','3',5],['2','3',10]]
    g1.instert_array_graph(array_g1,'0','3')


    first_level=graphs()
    first_arr=[['S', '1', 18], ['S', '2', 5], ['3', '4', 12], ['S', '3', 3], ['1', '6', 20], ['1', '7', 10], ['2', '5', 20], ['4', '8', 30], ['5', '9', 5], ['5', '10', 4], ['6', '10', 12], ['7', '11', 11], ['7', '12', 13], ['8', '11', 14], ['9', '10', 6], ['10', '13', 7], ['10', '14', 14], ['11', '15', 9], ['11', '12', 8], ['12', '16', 15], ['13', 'T', 2], ['14', 'T', 19], ['15', '16', 11], ['16', 'T', 8]]
    first_level.instert_array_graph(first_arr,'S','T')



    second_arr=[['T', '0', 5], ['T', '1', 3], ['T', '2', 4], ['T', '3', 6], ['T', '4', 7], ['T', '5', 4], ['T', '6', 1], ['T', '7', 6], ['T', '8', 7], ['T', '9', 5], ['36', 'S', 18], ['17', 'S', 6], ['60', 'S', 14], ['98', 'S', 11], ['49', 'S', 14], ['52', 'S', 9], ['56', 'S', 19], ['90', 'S', 16], ['15', 'S', 17], ['37', 'S', 15], ['0', 'S', 8], ['0', '91', 4], ['0', '72', 9], ['0', '50', 11], ['0', '20', 13], ['0', '75', 6], ['0', '47', 9], ['0', '86', 4], ['0', '44', 13], ['0', '16', 15], ['0', '54', 3], ['0', '38', 18], ['0', '11', 19], ['0', '19', 10], ['0', '65', 2], ['0', '37', 13], ['0', '67', 3], ['0', '36', 12], ['0', '61', 12], ['0', '31', 13], ['0', '48', 1], ['0', '77', 14], ['0', '70', 13], ['0', '64', 15], ['0', '17', 1], ['0', '92', 15], ['0', '89', 3], ['0', '27', 7], ['0', '55', 2], ['0', '85', 2], ['32', 'S', 5], ['28', 'S', 11], ['12', 'S', 18], ['93', 'S', 8], ['24', 'S', 15], ['74', 'S', 16], ['33', 'S', 4], ['75', 'S', 2], ['54', 'S', 19], ['21', 'S', 11], ['67', 'S', 19], ['5', 'S', 16], ['5', '14', 4], ['5', '71', 4], ['71', 'S', 3], ['81', 'S', 11], ['96', 'S', 14], ['59', 'S', 2], ['77', 'S', 16], ['27', 'S', 3], ['20', 'S', 10], ['41', 'S', 11], ['66', 'S', 14], ['88', 'S', 16], ['26', 'S', 16], ['14', 'S', 7], ['97', 'S', 18], ['47', 'S', 18], ['22', 'S', 9], ['38', 'S', 13], ['34', 'S', 11], ['79', 'S', 9], ['7', 'S', 1], ['7', '59', 3], ['7', '10', 3], ['7', '83', 17], ['7', '97', 3], ['7', '51', 15], ['50', 'S', 4], ['99', 'S', 11], ['99', '16', 16], ['99', '19', 1], ['99', '33', 16], ['99', '48', 11], ['99', '29', 1], ['99', '58', 3], ['99', '7', 4], ['99', '17', 10], ['99', '11', 6], ['99', '62', 11], ['99', '1', 19], ['99', '44', 16], ['99', '28', 8], ['99', '73', 1], ['99', '13', 14], ['99', '0', 7], ['99', '41', 7], ['99', '64', 12], ['99', '67', 11], ['99', '52', 18], ['99', '54', 5], ['99', '2', 3], ['99', '75', 10], ['99', '50', 6], ['99', '98', 8], ['99', '4', 7], ['99', '72', 3], ['99', '51', 5], ['99', '45', 16], ['99', '21', 12], ['99', '25', 16], ['99', '26', 17], ['99', '86', 3], ['99', '79', 9], ['99', '71', 9], ['99', '63', 6], ['99', '35', 3], ['99', '34', 9], ['99', '90', 16], ['99', '92', 10], ['99', '6', 8], ['99', '82', 14], ['99', '46', 17], ['99', '40', 11], ['99', '31', 11], ['99', '94', 5], ['99', '57', 12], ['99', '96', 16], ['99', '20', 17], ['99', '23', 14], ['99', '91', 11], ['99', '66', 7], ['99', '38', 12], ['51', 'S', 7], ['39', 'S', 16], ['55', 'S', 15], ['85', 'S', 9], ['72', 'S', 17], ['57', 'S', 3], ['63', 'S', 19], ['64', 'S', 1], ['6', 'S', 1], ['6', '82', 16], ['6', '22', 16], ['16', 'S', 2], ['16', '46', 10], ['1', '90', 17], ['1', '45', 11], ['1', '76', 1], ['1', '18', 7], ['1', '58', 12], ['1', '40', 11], ['2', '49', 12], ['2', '30', 8], ['2', '81', 6], ['2', '84', 7], ['2', '13', 6], ['2', '66', 1], ['2', '94', 15], ['2', '93', 11], ['2', '88', 6], ['2', '39', 19], ['2', '87', 8], ['3', '34', 13], ['3', '33', 10], ['3', '78', 3], ['3', '74', 6], ['3', '53', 11], ['3', '56', 4], ['3', '80', 6], ['3', '43', 10], ['3', '60', 3], ['3', '25', 12], ['3', '63', 15], ['3', '24', 4], ['4', '15', 11], ['4', '21', 19], ['4', '42', 9], ['4', '79', 14], ['4', '29', 6], ['4', '12', 3], ['4', '68', 3], ['4', '41', 10], ['4', '62', 10], ['4', '28', 7], ['4', '98', 11], ['4', '69', 19], ['4', '57', 12], ['4', '23', 9], ['4', '26', 12], ['8', '35', 14], ['8', '52', 8], ['9', '96', 5], ['9', '73', 12], ['9', '95', 3], ['9', '32', 11]]

    third_arr=[['T', '0', 5], ['T', '1', 7], ['T', '2', 1], ['T', '3', 3], ['T', '4', 16], ['T', '5', 7], ['T', '6', 3], ['T', '7', 4], ['T', '8', 4], ['T', '9', 5], ['T', '10', 3], ['T', '11', 3], ['126', 'S', 19], ['69', 'S', 2], ['45', 'S', 2], ['102', 'S', 1], ['63', 'S', 13], ['84', 'S', 9], ['24', 'S', 10], ['97', 'S', 6], ['138', 'S', 16], ['51', 'S', 19], ['65', 'S', 13], ['38', 'S', 13], ['52', 'S', 11], ['49', 'S', 9], ['95', 'S', 6], ['11', 'S', 7], ['11', '118', 19], ['11', '29', 11], ['80', 'S', 17], ['91', 'S', 15], ['115', 'S', 6], ['33', 'S', 1], ['67', 'S', 7], ['53', 'S', 13], ['98', 'S', 3], ['14', 'S', 12], ['127', 'S', 2], ['7', 'S', 6], ['7', '59', 12], ['62', 'S', 6], ['42', 'S', 11], ['140', 'S', 18], ['105', 'S', 10], ['99', 'S', 11], ['121', 'S', 18], ['54', 'S', 16], ['19', 'S', 16], ['120', 'S', 19], ['18', 'S', 7], ['25', 'S', 9], ['8', 'S', 16], ['124', 'S', 12], ['76', 'S', 18], ['46', 'S', 13], ['32', 'S', 11], ['82', 'S', 9], ['129', 'S', 14], ['29', 'S', 1], ['123', 'S', 10], ['40', 'S', 9], ['59', 'S', 19], ['47', 'S', 7], ['5', 'S', 19], ['5', '44', 19], ['5', '35', 15], ['5', '22', 9], ['5', '67', 15], ['5', '134', 15], ['5', '112', 17], ['5', '124', 17], ['5', '131', 5], ['5', '74', 3], ['5', '123', 6], ['5', '128', 5], ['5', '87', 7], ['72', 'S', 3], ['27', 'S', 5], ['90', 'S', 11], ['135', 'S', 5], ['16', 'S', 12], ['66', 'S', 14], ['34', 'S', 12], ['86', 'S', 13], ['113', 'S', 9], ['85', 'S', 3], ['87', 'S', 17], ['20', 'S', 6], ['93', 'S', 3], ['139', 'S', 3], ['10', 'S', 8], ['137', 'S', 3], ['117', 'S', 1], ['21', 'S', 3], ['0', '125', 9], ['0', '63', 2], ['0', '43', 7], ['0', '93', 14], ['0', '141', 11], ['0', '68', 11], ['0', '109', 17], ['0', '71', 19], ['0', '99', 16], ['0', '140', 3], ['0', '16', 19], ['0', '97', 11], ['0', '100', 19], ['0', '91', 13], ['0', '23', 17], ['0', '129', 3], ['0', '132', 13], ['0', '138', 10], ['0', '130', 5], ['0', '30', 17], ['0', '83', 3], ['0', '76', 8], ['0', '61', 16], ['0', '121', 2], ['0', '21', 3], ['0', '79', 15], ['0', '25', 15], ['0', '52', 6], ['0', '113', 8], ['0', '86', 13], ['0', '88', 6], ['0', '90', 16], ['0', '28', 16], ['0', '20', 18], ['0', '64', 19], ['0', '101', 7], ['0', '53', 18], ['0', '36', 17], ['0', '15', 17], ['0', '46', 17], ['0', '57', 8], ['0', '17', 14], ['0', '66', 6], ['0', '47', 2], ['0', '58', 10], ['0', '116', 13], ['0', '127', 11], ['0', '114', 4], ['0', '14', 7], ['1', '33', 4], ['1', '38', 5], ['1', '31', 8], ['1', '26', 19], ['1', '18', 10], ['1', '78', 3], ['1', '102', 2], ['1', '42', 10], ['1', '82', 5], ['1', '126', 9], ['1', '135', 6], ['1', '72', 13], ['1', '41', 12], ['1', '95', 14], ['1', '45', 1], ['1', '80', 12], ['1', '12', 17], ['1', '65', 8], ['2', '24', 8], ['3', '73', 4], ['3', '105', 16], ['3', '62', 9], ['3', '85', 1], ['3', '40', 2], ['3', '107', 17], ['3', '27', 8], ['3', '34', 5], ['3', '119', 4], ['3', '69', 16], ['3', '37', 10], ['3', '98', 18], ['3', '115', 12], ['3', '39', 14], ['3', '136', 11], ['3', '55', 9], ['3', '51', 9], ['3', '89', 11], ['3', '103', 9], ['3', '122', 5], ['3', '92', 5], ['4', '96', 14], ['4', '50', 4], ['4', '75', 14], ['4', '120', 17], ['4', '111', 7], ['4', '70', 2], ['4', '56', 14], ['4', '117', 3], ['4', '19', 17], ['4', '49', 16], ['4', '94', 14], ['4', '81', 2], ['4', '84', 2], ['4', '32', 7], ['4', '110', 18], ['4', '13', 10], ['4', '106', 12], ['6', '77', 11], ['6', '108', 14], ['6', '54', 18], ['6', '60', 5], ['6', '48', 8], ['6', '137', 11], ['9', '104', 16], ['9', '142', 3], ['9', '139', 19], ['12', '133', 6], ['143', '105', 17], ['143', '118', 7], ['143', '44', 4], ['143', '3', 1], ['143', '120', 12], ['143', '98', 19], ['143', '1', 13], ['143', '70', 7], ['143', '82', 15], ['143', '125', 7], ['143', '45', 13], ['143', '61', 17], ['143', '111', 13], ['143', '100', 10], ['143', '20', 9], ['143', '51', 7], ['143', '141', 6], ['143', '131', 2], ['143', '27', 11], ['143', '46', 14], ['143', '80', 6], ['143', '14', 1], ['143', '130', 2], ['143', '38', 14], ['143', '79', 2], ['143', '30', 15], ['143', '117', 13], ['143', '4', 12], ['143', '13', 4], ['143', '26', 17], ['143', '10', 19], ['143', '24', 5], ['143', '19', 8], ['143', '110', 2], ['143', '115', 8], ['143', '81', 10], ['143', '2', 14], ['143', '59', 18], ['143', '90', 3], ['143', '50', 12], ['143', '73', 13], ['143', '65', 1], ['143', '113', 13], ['143', '75', 15], ['143', '76', 4], ['143', '107', 16], ['143', '85', 12], ['143', '11', 16], ['143', '56', 16], ['143', '87', 12], ['143', '133', 8], ['143', '112', 13], ['143', '136', 11], ['143', '6', 10], ['143', '129', 7], ['143', '135', 5], ['143', '25', 2], ['143', '5', 13], ['143', '23', 2], ['143', '40', 10], ['143', '95', 7], ['143', '28', 5], ['143', '55', 7], ['143', '77', 10], ['143', '89', 12], ['143', '34', 18], ['143', '128', 13], ['143', '109', 16], ['143', '52', 11], ['143', '54', 8], ['143', '15', 19], ['143', '37', 18], ['143', '64', 11], ['143', '140', 2], ['143', '108', 9], ['143', '92', 3], ['143', '97', 4], ['143', '48', 3]]

    print('Πρώτο επίπεδο \n')
    all_methods(first_arr,'S','T')
    print('------------------\n')
    print('Δεύτερο επίπεδο \n')
    all_methods(second_arr,'T','S')
    print('------------------\n')
    print('Τρίτο επίπεδο \n')
    all_methods(third_arr,'T','S')
                
        
        
