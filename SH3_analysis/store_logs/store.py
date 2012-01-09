import csv
import sqlite3
import time
import glob

def store_logs()
    t = time.time()

    # establish db connection
    conn = sqlite3.connect('sh3_logs.db')
    
    # get cursor to write to the db
    c = conn.cursor()

    # create the table if it does not exist in the database
    c.execute('''create table if not exists sh3 (replica_num integer, sequence_num integer, w real, w_nominal integer)''')
    
    for csvfile in glob.glob("*.csv"):
        csvdata = csv.reader(open(csvfile, "rb"))
        print "loading", csvfile
        
        for row in csvdata:
            # skip the comment header
            if csvdata.line_num == 1:
                continue
            c.execute('INSERT OR IGNORE INTO sh3 VALUES (?,?,?,?)', row)
        
        conn.commit()
    
        print "\n Time Taken: %.3f sec" % (time.time()-t)
    
    conn.commit()
    c.close()

def make_xtc_name(replica_num, sequence_num):
    return 'lgw' + str(replica_num) + str(sequence_num) + '.xtc'
    
def construct_filenames():
    for csvfile in glob.glob("*.csv"):
        csvdata = csv.reader(open(csvfile, "rb"))
        print "loading", csvfile
        
        for (replica_num, sequence_num, w, w_nominal) in csvdata:
            # skip the comment header
            if csvdata.line_num == 1:
                continue
            
            # this is beta = 1/(k*T), k is the boltzman constant
            beta = float(w)
            if beta > 1.67 and beta < 1.68:
                fileslist.append(make_xtc_name(replica_num, sequence_num))

     
if __name__ == '__main__':
    main()



