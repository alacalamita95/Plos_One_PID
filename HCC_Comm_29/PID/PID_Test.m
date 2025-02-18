function PID = PID_Test(X, C_tot)

    % inizialize JMI, RED, SYN
    [ncombtot, ~] = size(C_tot);
    % names of the table
    names = {'JMI','MI1','MI2','RED','U1','U2','SYN','II'};
    PID = zeros(ncombtot, numel(names));
    
    JMI = zeros(ncombtot, 1);
    MI1 = zeros(ncombtot, 1);
    MI2 = zeros(ncombtot, 1);
    RED = zeros(ncombtot, 1);
    U1 = zeros(ncombtot, 1);
    U2 = zeros(ncombtot, 1);
    SYN = zeros(ncombtot, 1);
    II = zeros(ncombtot, 1);

    parfor i = 1:ncombtot
        t = X(:, C_tot(i, 1));
        d1 = X(:, C_tot(i, 2));
        d2 = X(:, C_tot(i, 3));
    
        % Joint Mutual Information between the target t and 2 drivers d1,d2
        JMI(i) = gcmi_cc(t,[d1 d2]);
        
        % Mutual Information seperately between t and d1 and d2
        MI1(i) = gcmi_cc(t,d1);
        MI2(i) = gcmi_cc(t,d2);
        
        % Define the Redundancy as minumum of the two MI terms
        RED(i) = min(MI1(i),MI2(i));
        
        % Unique information carried out by the single
        U1(i) = MI1(i)-RED(i);
        U2(i) = MI2(i)-RED(i);
        
        % defining Synergy as the information left of the JMI if we remove the RED
        % term and seperately the information carried out by the single terms U1 U2
        SYN(i) = JMI(i)-RED(i)-U1(i)-U2(i);
        
        % interaction TE is actually a measure of the 
        % ‘net’ synergy manifested in the transfer of information from the two sources to the target.
        II(i) = JMI(i)-MI1(i)-MI2(i);
    end

    PID(:, 1) = JMI;
    PID(:, 2) = MI1;
    PID(:, 3) = MI2;
    PID(:, 4) = RED;
    PID(:, 5) = U1;
    PID(:, 6) = U2;
    PID(:, 7) = SYN;
    PID(:, 8) = II;
    PID = array2table(PID, 'VariableNames', names);

end