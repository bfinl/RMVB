function rad1 = find_ellipradius(A, dir)
dir = dir / norms(dir);
[U, S, ~]  = svds(A, rank(A));
rad1 = 1 / norms((S^-1) * U.' * dir);
end