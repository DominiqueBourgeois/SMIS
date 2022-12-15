function N = check_in_fov(x,y,f,im_par)
%Check whether molecules are still within FOV
%Use a margin factor of f
n_max=f*im_par.n;
m_max=f*im_par.m;
n_min=-(f-1)*im_par.n;
m_min=-(f-1)*im_par.m;

w=find(x>=n_min & y>=m_min & x<=n_max & y<=m_max);
N=length(w);

