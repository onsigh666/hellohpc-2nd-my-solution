## HelloHPC 2nd Writeup 

### 2.小交问答  
1. B. LibreChat 不是  

2. 4.3 3.4
3. 838 TFLOPS 104.8TFLOPS  use nvidia-smi or ncu to analyse.
4. 
    2: **enp189s0f0** : <BROADCAST,MULTICAST,UP,LOWER_UP> mtu 1500 qdisc fq_codel state UP mode DEFAULT group default qlen 1000  
    link/ether 10:c3:ab:50:cb:28 brd ff:ff:ff:ff:ff:ff  
    3: **enp189s0f1**: <NO-CARRIER,BROADCAST,MULTICAST,UP> mtu 1500 qdisc fq_codel state DOWN mode DEFAULT group default qlen 1000
    link/ether 10:c3:ab:50:cb:29 brd ff:ff:ff:ff:ff:ff  
    4: **enp125s0f0**: <NO-CARRIER,BROADCAST,MULTICAST,UP> mtu 1500 qdisc mq state DOWN mode DEFAULT group default qlen 1000
    link/ether 10:c3:ab:50:cb:27 brd ff:ff:ff:ff:ff:ff  
    **firmware-version: 1.8.15.0**  
      
    ***How to find the info:***
    First allocate a calculation node: `salloc -p kp_run -q kp_run -c 1 -t 00:05:00`  
    Then list all ip link: `ip link show`  
    Finally find detail information: `ethtool -i enp189s0f0`  

5.   Entry point address:               0x400580  
    In calculation node: `gcc --version`  
    Write a simple c code, like hello.c, and compile it.  
    `gcc hello.c -o hello`  
    `readelf -h hello | grep "Entry point address"`  
    The output is the entrance.

7.  f4e9f91cc7ffb0e4810744776728a8d863623bc518ab754ce90b6d1513720494
    ```
    [stu1562@kp032 ~]$ wget https://xflops.sjtu.edu.cn/xflops-htop.png -O xflops-htop.png  
    --2026-09-23 19:50:43--  https://xflops.sjtu.edu.cn/xflops-htop.png
    Resolving proxy1.pi.sjtu.edu.cn (proxy1.pi.sjtu.edu.cn)... 172.16.0.103
    Connecting to proxy1.pi.sjtu.edu.cn (proxy1.pi.sjtu.edu.cn)|172.16.0.103|:3128... connected.
    Proxy request sent, awaiting response... 200 OK
    Length: 460 [text/html]
    Saving to: ‘xflops-htop.png’

    xflops-htop.png 100%[====>]     460  --.-KB/s    in 0s      

    2026-09-23 19:50:43 (18.2 MB/s) - ‘xflops-htop.png’ saved [460/460]

    [stu1562@kp032 ~]$ ls
    hello  hello.c  main.c  test.slurm  xflops-htop.png
    [stu1562@kp032 ~]$ sha256sum xflops-htop.png 
    f4e9f91cc7ffb0e4810744776728a8d863623bc518ab754ce90b6d1513720494  xflops-htop.png
    ```
