History
========

When I was doing a series of Hi-C projects, I tried many drawing software such as 
HiTC, HiCPlotter, etc., and used python's matplotlib package to practice various Hi-C 
and multi-omics visualization drawing methods. As early as 2018, when I was doing a project about 
the relationship between pseudogene and three-dimensional genome, I had the idea of building a 
flexible multi-omics visualization software. And sorted out the original Hi-C triangular rectangle heatmap 
code at that time. In subsequent projects, we have been perfecting the implementation methods of 
various multi-omics tracks. For example, the Beizel curve is used to replace the ellipse to realize 
the drawing method of the loop, so that the connection of the loop looks smoother. Use the 
selection of X, Y, Z coordinates in pcolormesh to realize the triangular heat map to replace the original 
layered implementation method from the Hi-C mode. The triangular heatmap implementation method from HiTC was first 
converted to the python version by Lixiang Ma in our team. 


which project uses trackc
----------------

- PIBF1
- pseudogene
- zhongshan eye
- muxu
