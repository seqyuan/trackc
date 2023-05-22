History
========

While working on a series of Hi-C projects, I experimented with various Hi-C data drawing software such as 
HiTC and HiCPlotter, as well as utilized Python's Matplotlib package to practice different Hi-C and 
multi-omics visualization techniques. As early as 2018, during a project exploring the relationship 
between pseudogenes and 3D genome, I had an idea to create a flexible multi-omics visualization software. 
I began by sorting out the original Hi-C triangular rectangle heatmap code and continued refining its 
implementation in subsequent projects to include various multi-omics tracks. For instance, we replaced the 
ellipse with Bezier curves to create smoother loop connections, and used X, Y, and Z coordinate selection in 
pcolormesh to implement the triangular heatmap instead of the original layered method from the Hi-C mode. 
The implementation method for the triangular heatmap from HiTC was initially converted to Python by Lixiang Ma 
one our team member.

which project uses trackc
-------------------------

- PIBF1
- pseudogene
- zhongshan eye
- muxu
