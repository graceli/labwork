import ssh

class Cluster():
    """
    Cluster
    
    This is a compute resource that can be accessed via SSH. 
    This class is the base class for which all clusters inerhit their function.
    """

    def _executeCommand(self, remote_command):

        client = ssh.SSHClient()
        client.load_system_host_keys()
        client.connect(self.host, username=self.user)

        #this returns a tuple of stdin, stdout, stderr
        stdin, stdout, stderr = returned_data = client.exec_command(remote_command)
        #stdin, stdout, stderr = returned_data = client.exec_command("source ~/.bash_profile; " + remote_command)
        exit_status = stdout.channel.recv_exit_status()

        return stderr.readlines(), stdout.readlines()

class Orca(Cluster):
    """
    Orca Cluster

    This is a subcluster of Sharcnet
    """

    def __init__(self, user):
        self.user = user
        self.name = "Orca"
        self.host = "orca.sharcnet.ca"
        self.submit_command = "sqsub"
        #self.joblist_command = "/opt/sharcnet/sq-tm/2.4/bin/sqjobs | grep ceing"
        self.joblist_command = "/opt/sharcnet/util/current/bin/gjobs | grep orca"

if __name__ == '__main__':
    x = Orca("ceing")
    print x._executeCommand(x.joblist_command)
    #print x._executeCommand("du -h")

