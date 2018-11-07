import socket
import os

server_adress = "./udp_socket"
sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)

sock.bind()
sock.listen(1)
connection, client_adress = sock.accept()
message = connection.recv(1024)
print(message)
connection.close()

