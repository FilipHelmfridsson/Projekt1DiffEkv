function resultatEuler = euler(u, N, tslut)

    h = tslut / N;

    for i = 1:N
        u = [u, u(:,end) + h * planetrhs(u(:,end))];
    end

    resultatEuler = u';

end