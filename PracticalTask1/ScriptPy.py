
import shlex, subprocess

log = open('rep.txt','w')

for i in range(6):
    res = 0.0
    for j in range(10):
        c = subprocess.check_output(["./main","_rA","_rB","_rC",str(i)])
        l = c.decode("utf-8")
	res = res + int(l)
    res = res / 10
    log.write(str(i))
    log.write(" ")
    log.write(str(res)+'\n')
    print (res)
subprocess.Popen("./gnuScript.sh")