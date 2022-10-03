# Purpose: Place spray gun(s) over discretized square grid to get
# composition profile of an experiment

# Credit to Steven Torrisi's CompositionSpray code for inspiration!

import scipy
from pymatgen.core import Composition, Element
from typing import List
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

class ElementSource(object):
    
    def __init__(self,x0:float,y0:float,sigma:float,element:Element, intensity: float = 1.):
        '''
        Args:
            x0, y0 : Coordinates of spray gun
            sigma : Standard Deviation of Gaussian Spray Profile (NOTE: height dependent?)
            element : Element of the PeriodicTable associated with the spray gun
            intensity : Effectively concentration of spray
        '''
        self.x0 = x0
        self.y0 = y0
        self.sigma = sigma
        self.element = element
        self.intensity = intensity
        

    def pdf(self,x,y):
        
        #prefactor = 1/(self.sigma **2 * 2 * np.pi)
        prefactor = 1
        return self.intensity * prefactor \
            * np.exp( -(x-self.x0)**2 / (2* self.sigma **2))\
            * np.exp( -(y - self.y0)**2 /(2*self.sigma **2))

class GridExperiment(object):

    def __init__(self,x:float,y:float,\
    x_steps:int=10,y_steps:int=10, spray_guns:List[ElementSource] = None):
        self.x = x
        self.y = y
        self.x_steps = x_steps
        self.y_steps = y_steps
        self.spray_guns = spray_guns if spray_guns else []

        self.df = pd.DataFrame()

    def comp_to_color(self,comp:Composition,comp_to_color:dict = None)->tuple:
    
    
        compdict = comp.as_dict()
        elts = sorted(list(compdict.keys()))
        #print(elts)
        values = [compdict[key] for key in elts]
        norm = sum(values)
    
        color_tup = tuple(values/norm)

        return color_tup

    def combine_element_sources(self,x:float,y:float, sources: List[ElementSource],normalize: bool = True)->Composition:
    
        elements = {}
        
        #print(sources)
        for src in sources:
            if src.element in elements:
                elements[src.element] += src.pdf(x,y)
            else:
                elements[src.element] = src.pdf(x,y)
         
        if normalize:
            elt_sum = np.sum(list(elements.values()))
            
            for key in elements.keys():
                elements[key] /= elt_sum
        
        return Composition(elements)

    def run(self,plot_graph=True):
        '''
        Run the experiment, return array of data points, plot graph if True
        '''
        #element_list = np.unique([str(sg.element) for sg in self.spray_guns])
        points = []
        color = []
        for x in np.linspace(0,self.x,self.x_steps):
            for y in np.linspace(0,self.y,self.y_steps):
                comp = self.combine_element_sources(x,y,self.spray_guns,normalize=True)
                color.append(self.comp_to_color(comp))
                points.append({'x':x,'y':y} | comp.as_dict())
        
        df = pd.DataFrame(points)
        df = df.replace(np.nan,0)

        color_arr = [f'rgba(0,{x[0]},{x[1]},.7)' for x in color] 
        if plot_graph:
            #alpha = .6
            #colors = [
            #    [136, 0, 255, alpha],
            #    [255, 243, 138, alpha],
            #    [143, 251, 255, alpha],
            #    [255, 161, 179, alpha],
            #]

            size_df = df.loc[:,~df.columns.isin(["x","y"])].sum(axis=1)
            max_size = max(size_df)
            size_df=size_df/max_size

            fig = px.scatter(df,x="x",y="y",hover_data=df.columns,width=800,height=800)
            fig.update_traces(marker={'size': 20, 'color':color_arr})
            fig.update_layout(
                #title="Plot Title",

                xaxis_title="X ($X (\mu m)$)m)",
                yaxis_title="Y ($Y (\mu m)$",

                #legend_title="Legend Title",
                font=dict(
                    family="Arial, monospace",
                    size=24,
                    #color="RebeccaPurple"
                )
            )
            fig.show()
        
        return df
                

