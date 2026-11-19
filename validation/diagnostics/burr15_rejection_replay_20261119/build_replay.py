"""Offline snapshot replay linked to unchanged frozen 8a85e37 objects."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess

repo = Path('/home/leonardoab/Documentos/PRIMS-F-burr15-20261119')
frozen = repo.parent / 'PRIMS-F-burr15-regression-8a85e37-20261119'
runs = repo.parent / 'PRIMS-F-burr15-results-20261119/v1/runs'
report = repo / 'validation/diagnostics/burr15_rejection_replay_20261119'
report.mkdir(parents=True, exist_ok=True)
build = report / 'build'
build.mkdir(exist_ok=True)
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
source = (repo / 'validation/diagnostics/bg_ss1_stall/probe_bg_ss1_stall.f90').read_text()
replacements = [
    ('integer,parameter::nc=52,np=48', 'integer,parameter::nc=106,np=102'),
    ('norm_saved,theta,width,dummy(3)', 'norm_saved,theta,width,dummy(3),rsaved(3,np)'),
    ('candidate(nc,3),lambda', 'candidate(nc,3),lambda,lambda0,nt,ps_si,jl_si,jg_si,mul_si,mug_si'),
    ('info,pointfails', 'info,pointfails,firstbad,firststatus,jacmode,firstface,firstfacestatus'),
    ('  cfg%omega_u=.53_wp;cfg%omega_rho=997.05_wp;cfg%omega_P=101325.0_wp',
     "  open(newunit=iu,file='probe_parameters.dat',status='old')\n"
     '  read(iu,*)cfg%omega_u,ps_si,jl_si,jg_si,mul_si,mug_si,jacmode\n'
     '  close(iu)\n'
     '  cfg%omega_rho=997.05_wp;cfg%omega_P=101325.0_wp'),
    ('208000.0_wp/cfg%omega_P', 'ps_si/cfg%omega_P'),
    ('528.4365_wp/(cfg%omega_rho*cfg%omega_u)', 'mul_si/(cfg%omega_rho*cfg%omega_u)'),
    ('.13244340430539120_wp/(cfg%omega_rho*cfg%omega_u)', 'mug_si/(cfg%omega_rho*cfg%omega_u)'),
    ('.53_wp/cfg%omega_u;cfg%njg=.11_wp/cfg%omega_u', 'jl_si/cfg%omega_u;cfg%njg=jg_si/cfg%omega_u'),
    ('cfg%jac_dissipation_fd=.false.', 'cfg%jac_dissipation_fd=jacmode==1'),
    ("cfg%junction_model='MATCHED';", "cfg%gas_eos='BAROTROPIC';cfg%junction_model='MATCHED';"),
    ('cfg%N_pipe=24', 'cfg%N_pipe=51'),
    ('3,50', '3,104'), ('3,49', '3,103'),
    ('(51,:)', '(105,:)'), ('(52,:)', '(106,:)'), ('(50,:)', '(104,:)'),
    ('(51:52)', '(105:106)'), ('(50)', '(104)'),
    ('(26,:)', '(53,:)'), ('(27,:)', '(54,:)'), ('(26)', '(53)'), ('(27)', '(54)'),
    ('(:,:,26)', '(:,:,53)'),
    ('    cfg%thetav(i)=theta;cfg%cell_width(i)=width',
     '    cfg%thetav(i)=theta;cfg%cell_width(i)=width;rsaved(:,i-2)=dummy'),
    ('dissipation_matrix=sim%lagged_dissipation(:,:,i))',
     'dissipation_matrix=sim%lagged_dissipation(:,:,i),left_cell=i)'),
    ('dissipation_matrix=sim%lagged_dissipation(:,:,53))',
     'dissipation_matrix=sim%lagged_dissipation(:,:,53),left_cell=53)'),
    ("    write(*,'(a,i3,i8,es25.16)')'assembly Jacobian profile,status,norm ',profile,code,norm",
     "    write(*,'(a,i3,i8,es25.16)')'assembly Jacobian profile,status,norm ',profile,code,norm\n"
     "    if(profile==2.and.code==0)write(*,'(a,es25.16)')'residual difference from snapshot ',maxval(abs(r-rsaved))"),
    ('      do k=0,4\n        lambda=2.0_wp**(-k);candidate=state',
     '      lambda0=1\n'
     '      do i=3,104\n'
     '        do k=1,2\n'
     '          if(dx(k,i-2)<0)lambda0=min(lambda0,-0.90_wp*state(i,k)/dx(k,i-2))\n'
     '        enddo\n'
     '      enddo\n'
     '      do k=0,17\n        lambda=lambda0*2.0_wp**(-k);candidate=state'),
    ('        pointfails=0',
     '        call sim%assemble(candidate,old,dt,r,ab,bb,cb,nt,fin,fout,code,.false.)\n'
     "        write(*,'(a,i3,es25.16,i8,es25.16,l2)')'proposal k,lambda,status,norm,Armijo ',k,lambda,code,nt, &\n"
     '          code==0.and.(nt<norm*(1-1e-4_wp*lambda).or.nt<=cfg%tol_res)\n'
     '        pointfails=0;firstbad=0;firststatus=0'),
    ('          if(code/=0)pointfails=pointfails+1',
     '          if(code/=0)then\n'
     '            pointfails=pointfails+1\n'
     '            if(firstbad==0)then\n'
     '              firstbad=i;firststatus=code\n'
     '            endif\n'
     '          endif'),
    ("        call junction_probe(candidate,10+k)",
     "        write(*,'(a,2i8)')'first point failure cell,status ',firstbad,firststatus\n"
     '        firstface=0;firstfacestatus=0\n'
     '        if(pointfails==0)then\n'
     '          do i=3,103\n'
     '            call sim%face_pair(candidate(i,:),candidate(i+1,:),cfg%thetav(i),cfg%thetav(i+1),fm,fp,code, &\n'
     '              dissipation_matrix=sim%lagged_dissipation(:,:,i),left_cell=i)\n'
     '            if(code/=0)then\n'
     '              firstface=i;firstfacestatus=code;exit\n'
     '            endif\n'
     '          enddo\n'
     "          write(*,'(a,2i8)')'first face failure leftcell,status ',firstface,firstfacestatus\n"
     '        endif\n'
     '        call junction_probe(candidate,10+k)'),
    ('    call sim%point_state(uu(54,:),cfg%thetav(54),vr,fr,source,aa,ss,st)',
     "    if(st/=0)then\n      write(*,'(a,2i8)')'junction left cell status ',iprofile,st\n      return\n    endif\n"
     '    call sim%point_state(uu(54,:),cfg%thetav(54),vr,fr,source,aa,ss,st)\n'
     "    if(st/=0)then\n      write(*,'(a,2i8)')'junction right cell status ',iprofile,st\n      return\n    endif"),
]
for before, after in replacements:
    assert before in source, before
    source = source.replace(before, after)
probe = report / 'probe_rejection_v2.f90'
probe.write_text(source)
objects = sorted(path for path in (frozen / 'tmp/release').glob('*.o') if path.name != 'main.o')
command = [
    'gfortran', '-std=f2018', '-O2', '-ffree-line-length-none', '-fopenmp',
    '-I' + str(frozen / 'tmp/release'), '-J' + str(build), str(probe),
    *map(str, objects), '-llapack', '-lblas', '-o', str(build / 'probe_rejection'),
]
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
process = subprocess.run(command, capture_output=True, text=True, env=env)
(build / 'build.log').write_text(process.stdout + process.stderr)
assert process.returncode == 0, process.stderr
provenance = {
    'frozen_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=frozen, text=True).strip(),
    'source_worktree': str(frozen), 'command': command,
    'source_sha256': {str(p.relative_to(frozen)): sha(p) for p in (frozen / 'scr').glob('*.f90')},
    'objects_sha256': {p.name: sha(p) for p in objects},
    'probe_sha256': sha(probe), 'probe_binary_sha256': sha(build / 'probe_rejection'),
    'frozen_binary_sha256': sha(frozen / 'bin/prims_f'), 'runs': {},
}
for case in ('D_SS1_fixed', 'D_SS3_fixed'):
    directory = runs / case
    manifest = json.loads((directory / 'manifest.json').read_text())
    assert provenance['frozen_head'] == manifest['git_head']
    assert all(provenance['source_sha256'][name] == value for name, value in manifest['source_sha256'].items())
    controls = manifest['controls']
    scales = (directory / 'run_metadata.txt').read_text().split('Lr,CA,Cg,Gl,Gg,omega_P,omega_rho,omega_u:')[1].splitlines()[0].split()
    item = {'snapshot_sha256': sha(directory / 'last_rejected_state.dat'), 'original_manifest': manifest, 'replays': {}}
    for flag in (0, 1):
        target = report / (case + ('_F_v2' if flag == 0 else '_T_v2'))
        target.mkdir(exist_ok=True)
        shutil.copyfile(directory / 'last_rejected_state.dat', target / 'last_rejected_state.dat')
        parameters = [float(scales[-1]), controls['ps'], controls['jl'], controls['jg'], controls['mul'], controls['mug']]
        (target / 'probe_parameters.dat').write_text(' '.join(format(x, '.17e') for x in parameters) + f' {flag}\n')
        result = subprocess.run([str(build / 'probe_rejection')], cwd=target, env=env, capture_output=True, text=True, timeout=120)
        (target / 'probe.log').write_text(result.stdout)
        (target / 'stderr.log').write_text(result.stderr)
        assert result.returncode == 0, result.stderr
        item['replays'][str(flag)] = {'returncode': result.returncode, 'directory': str(target), 'log_sha256': sha(target / 'probe.log')}
        print(case, 'FT'[flag], result.returncode, flush=True)
    provenance['runs'][case] = item
(report / 'provenance_v2.json').write_text(json.dumps(provenance, indent=2) + '\n')
