import ftplib

def init_ftp(host, urs, pwd):
    ftp = ftplib.FTP(host, user, pwd)
    ftp.encoding = "utf-8"
    return ftp

def download(ftp, localpath, ftppath):
    file = open(localpath, "wb")
    ftp.retrbinary("RETR " + ftppath, file.write)
    file.close()