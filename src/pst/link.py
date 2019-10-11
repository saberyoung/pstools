"""############################################################################ 
2019/1/30 Start
define some functions for linking soft/email/etc
""" ############################################################################
from __future__ import print_function
from builtins import input
import os

def sendemail(email,emailpass,emailsmtp,subj,\
              fromaddr,toaddrs,text,Imglist=[],Filelist=[]):
    """
    send email
    """
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.image import MIMEImage
    from email.mime.multipart import MIMEMultipart, MIMEBase
    from email import encoders

    try:
        # settings
        msg = MIMEMultipart()
        msg['Subject'] = subj
        msg['From'] = fromaddr
        msg['To'] = toaddrs

        # text
        text = MIMEText(text)
        msg.attach(text)

        # attachments
        for FileName in Filelist:
            part = MIMEBase('application', "octet-stream")
            part.set_payload(open(FileName, "rb").read())
            encoders.encode_base64(part)
            part.add_header('Content-Disposition', 'attachment; filename="%s"'%FileName)
            msg.attach(part)

        # image
        for ImgFileName in Imglist:
            img_data = open(ImgFileName, 'rb').read()
            image = MIMEImage(img_data, name=os.path.basename(ImgFileName))
            msg.attach(image)

        server = smtplib.SMTP(emailsmtp)
        server.starttls()
        server.login(email,emailpass)
        server.sendmail(fromaddr,toaddrs,msg.as_string())
        server.quit()
        return True
    except:
        return False

def createSSHClient(server, port, user, password):

    import paramiko
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, user, password)
    return client

def sftp_put_dir(server, port, user, password, filelist, remote_dir):

    import paramiko    
    t = paramiko.Transport(sock=(server, port))
    t.connect(username=user, password=password)
    sftp = paramiko.SFTPClient.from_transport(t)   

    try:sftp.chdir(remote_dir)  # Test if remote_path exists
    except IOError:sftp.mkdir(remote_dir)  # Create remote_path       

    for x in filelist:
        filename = os.path.split(x)[-1]
        remote_filename = remote_dir + '/' + filename
        print('Put %s...' % filename)
        sftp.put(x, remote_filename)
    sftp.close()

def wechat(_name,_msg,_img):
    # pip install -U wxpy
    # pip install -U wxpy -i "https://pypi.doubanio.com/simple/"
    from wxpy import Bot

    # initialize
    bot = Bot()

    # find friend
    if len(bot.friends().search(_name)) == 1:answ = 0
    else:
        for _nn,_friend in enumerate(bot.friends().search(_name)):print(_nn+1,_friend)
        answ = int(input('which one?'))
    my_friend = bot.friends().search(_name)[answ-1]
    
    # send msg
    if _msg: my_friend.send(_msg)
    # send img
    if _img: my_friend.send_image(_img)

def slack(token, usr, content):
    # pip install slackclient   
    try: from slackclient import SlackClient
    except: return False

    slack_client = SlackClient(token)    
    slack_client.api_call("chat.postMessage", channel=usr,
                          text=content, as_user=True)
    return True

def phone(_account,_token,_from,_to,_txt):
    # pip install twilio
    try:from twilio.rest import Client
    except: return False

    account_sid = _account
    auth_token  = _token
    client = Client(account_sid, auth_token)
    message = client.messages.create(
        to=_to,
        from_=_from,
        body=_txt)
    print(message.sid)
    return True
