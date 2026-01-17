% replicate_original_results.m
function replicate_original_results()
    try
        cd('/mnt/data/replication/source/Code/Code for JME');
        rng(123);
        if ~exist('/mnt/data/replication/outputs/logs', 'dir'), mkdir('/mnt/data/replication/outputs/logs'); end
        diary('/mnt/data/replication/outputs/logs/run_20250921_214513.txt');
        diary on;
        fprintf('\n=== Replication Run Started: %s ===\n', datestr(now));
        fprintf('CWD: %s\n', pwd);
        fprintf('Seed: 123\n');
        % See comments for centralizing fsolve options.
        run('main.m');
        outdir = '/mnt/data/replication/outputs';
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        save(fullfile(outdir, 'run_workspace.mat'));
        try
            if exist('Pc','var') && exist('C','var') && exist('Px','var') && exist('X','var') ...
               && exist('NX','var') && exist('w','var') && exist('L','var') ...
               && exist('r','var') && exist('K','var')
                bc_resid = max(abs(Pc.*C + Px.*X + NX - (w.*L + r.*K)), [], 'all');
                fprintf('\nMax budget-constraint residual: %g\n', bc_resid);
            else
                fprintf('\n[Note] Budget variables not all found; skipping BC residual check.\n');
            end
        catch ME
            fprintf('\n[Warning] Post-run BC check failed: %s\n', ME.message);
        end
        diary off;
        fprintf('\n=== Replication Run Finished ===\n');
    catch ME
        try, diary off; end
        rethrow(ME);
    end
end
