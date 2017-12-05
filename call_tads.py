#!/usr/bin/env python
import subprocess
import signal

# Call tads for 20 encode cell lines using matrix2insulation.pl and insulation2tads.pl

def get_chrom_matrices(chrom):
	s = subprocess.Popen('ls /data/projects/hic/TAD_consensus/borrman/*/*/C-40000/iced/*' + 
		chrom + '.matrix.gz', stdout=subprocess.PIPE, shell=True)
	cfiles = s.communicate()[0].strip().split('\n')
	return cfiles

def get_i_and_b(prefix):
	# Get insulation and boundary files
	s = subprocess.Popen('ls ' + prefix + '*.insulation', stdout=subprocess.PIPE, shell=True)
	s.wait()
	i = s.communicate()[0].strip()
	s = subprocess.Popen('ls ' + prefix + '*.boundaries', stdout= subprocess.PIPE, shell=True)
	s.wait()
	b = s.communicate()[0].strip()
	return i, b


def main():
	for c in map(str, range(1, 23)) + ['X']:
		# Insulation and boundary calls
		ps = []
		cm = get_chrom_matrices('chr' + c)
		for m in cm:
			m_prefix = m.strip().split('/')[-1][:-10]
			p = subprocess.Popen('perl ~/cworld-dekker/scripts/perl/matrix2insulation.pl -i ' +
				m + ' -v 2> ' + m_prefix + '_insulation.out', shell=True,
				 preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
			# signal module used to fix gzip: stdout: broken pipe bug as explained here:
			# https://blog.nelhage.com/2010/02/a-very-subtle-bug/
			ps.append(p)
		for p in ps:
			p.wait()
		
		#TAD calls
		ps2 = []
		for m in cm:
			m_prefix = m.strip().split('/')[-1][:-10]
			i, b = get_i_and_b(m_prefix)
			p2 = subprocess.Popen('perl ~/cworld-dekker/scripts/perl/insulation2tads.pl -i ' + i +
				' -b ' + b + ' -v 2> ' + m_prefix + '_tad.out', shell=True)
			ps2.append(p2)
		for p2 in ps2:
			p2.wait()


if __name__ == '__main__':
	main()