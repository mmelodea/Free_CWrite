------------------------
# Hello world
a = "Hello"
b = "world"
print a,b
------------------------

#Ler 2 listas e criar 3ª lista com elementos intercalados das 2 listas lidas
l=[1,2,3,4,5]
g=[6,7,8,9,10]
t=[]
for y in range(0,5):
   t.append(l[y])
   t.append(g[y])
   
print t
-----------------------------------------------------------


#Listas
meses = ['Janeiro','Fevereiro','Marco','Abril','Maio','Junho','Julho','Agosto','Setembro','Outubro','Novembro','Dezembro'] 
while 1:
    mes = input(“Escolha um mes (1-12)?”)
    if 1 <= mes <= 12:
print 'O mes é', meses[mes - 1]
----------------------------------------------------------


#Recebe 4 notas e retorna algumas informações
notas=[]
a=0

while 1:
    nota = input("Insira a nota: ")
    if nota != -1:
      notas.append(nota)
      a+=nota
    elif nota == -1:
      break

media=a/len(notas)
print"\n#####################################"
print "Quantidade de notas: %d" % len(notas)
print "Notas inseridas: ",notas
print "Media= %.2f" % (media)

upm=[]
for i in notas:
   if i > (media):
     upm.append(i)
     
print "Notas acima da media: %d" % len(upm)
------------------------------------------------------------

#Fatorial
resultado = 1
num = input("Entre com um numero inteiro: ")
num2 = num
while num2 > 1:
   resultado = resultado*num2
   num2 = num2 - 1
print num,"!: %.3e" % resultado
---------------------------------------------------------


#Funções
def juiz(reu):
     if reu > 0:
        print reu,": P"
     elif reu < 0:
        print reu," :N"
     else:
        print reu," : Indeterminavel!"
        
n = input("Insira um numero para julgamento: ")
juiz(n)
-------------------------------------------------------------


#Matricula identifica periodo e tempo restante de conclusão
def ID(m):
     l = []
     for i in m:
         l.append(i)
     
     if l[-2] == "1":
        print "Seu periodo: Primeiro"
        print "Tempo de formacao: ",8*4*30,"dias"
     elif l[-2] == "2":
        print "Seu periodo: Segundo"
        print "Tempo de formacao: ",7*4*30,"dias"
     elif l[-2] == "3":
        print "Seu periodo: Terceiro"
        print "Tempo de formacao: ",6*4*30,"dias"
     elif l[-2] == "4":
        print "Seu periodo: Quarto"
        print "Tempo de formacao: ",5*4*30,"dias"
     elif l[-2] == "5":
        print "Seu periodo: Quinto"
        print "Tempo de formacao: ",4*4*30,"dias"
     elif l[-2] == "6":
        print "Seu periodo: Sexto"
        print "Tempo de formacao: ",3*4*30,"dias"
     elif l[-2] == "7":
        print "Seu periodo: Setimo"
        print "Tempo de formacao: ",2*4*30,"dias"
     elif l[-2] == "8":
        print "Seu periodo: Oitavo"
        print "Tempo de formacao: ",1*4*30,"dias"

n = raw_input("Insira sua matricula entre '': ")
ID(n)
------------------------------------------------------------


#Recursividade de funções
>>>def fatorial(n):
if n <= 1:
...
return 1
...
return n * fatorial(n - 1)
...
...
>>>print ’2! = ’,fatorial(2)
2! = 2
-----------------------------------------------------------------


#Recebe valores do teclado e os dispõem em ordem crescente e decrescente
s=[]
f = input("Diga o tamanho da serie: ")
for i in range(0,f):
    e= input("Insira um elemento: ")
    s.append(e)

s.sort()    
print s

t=[]
for i in range(1,len(s)+1):
    t.append(s[-i])

print t
--------------------------------------------------------------


#Funções
def imprime_cardapio (lista_pratos):
  print "Cardapio para hoje\n"
  for p in lista_pratos:
    imprime_prato(p)
  print "\nTotal de pratos: %d" % len(lista_pratos)

def imprime_prato(p):
   print "%s ........ %10.2f" % (p["nome"], p["preco"])

p1 = {"nome": "Arroz com brocolis","preco": 9.90}
p2 = {"nome": "Soja com legumes","preco": 7.80}
p3 = {"nome": "Lentilhas", "preco": 4.80}
lista_pratos = [p1, p2, p3]

imprime_cardapio(lista_pratos)
----------------------------------------------------------

#abrir e manipular arquivo
a = open("arquivo.txt","r")
texto = a.read( )
print texto

saída: ... bla bla bla ...

a = open("arquivo.txt","w ou a")
a.write("bla bla bla")
---------------------------------------------------------

#entrada e saída de dados
a = raw_input("Entre com um numero de 0 a 10: ")
n = int(a)
if not 0<n<10:
 print "Numero invalido!"
elif n % 2 == 0:
 print "E par!"
else:
 print "E impar!"
---------------------------------------------------------

#construindo classes
class Retangulo:
       lado_a = None  #aparentemente, não usar isto, não gera erro no programa
       lado_b = None

def __init__(self, lado_a, lado_b):
       self.lado_a = lado_a
       self.lado_b = lado_b
       print "Criando nova instância Retangulo"

def calcula_area(self):
    return self.lado_a * self.lado_b

def calcula perímetro(self):
    return 2 * self.lado_a + 2 * self.lado_b
    
#Uso da classe acima:
# from file.py import Retangulo
# r=Retangulo(2,4)
# r.calcula_area( ) ou r.calcula_perimetro( )
-------------------------------------------------------------------


#Estado dos alunos de uma turma
alunos = ['Fred','Suzana','Claudio','Puga','Robson','Gustavo']
nota = [5.4, 6.2, 2.9, 9.9, 7.8, 4.9]
faltas = [9, 5, 15, 2, 11, 12]

contador = 0
for aluno in alunos:

 if nota[contador] >= 6.0 and faltas[contador] <= 10:
  print '\nAluno: ',aluno
  print 'Nota final: ',nota[contador]
  print 'Faltas: ',faltas[contador]
  print 'Resultado: Passou de ano'

 elif nota[contador] >= 6.0 and faltas[contador] > 10:
  print '\nAluno: ',aluno
  print 'Nota final: ',nota[contador]
  print 'Faltas: ',faltas[contador]
  print 'Resultado: Recuperacao por falta'

 elif nota[contador] >= 4.0 and nota[contador] < 6.0 and faltas[contador] <= 10:
  print '\nAluno: ',aluno
  print 'Nota final: ',nota[contador]
  print 'Faltas: ',faltas[contador]
  print 'Resultado: Recuperacao por nota'

 elif nota[contador] >= 4.0 and nota[contador] < 6.0 and faltas[contador] > 10:
  print '\nAluno: ',aluno
  print 'Nota final: ',nota[contador]
  print 'Faltas: ',faltas[contador]
  print 'Resultado: Repetiu direto por nao obter nota e por excesso de faltas'

 elif nota[contador] < 4.0:
  print '\nAluno: ',aluno
  print 'Nota final: ',nota[contador]
  print 'Faltas: ',faltas[contador]
  print 'Resultado: Repetiu direto por nota'
  
 contador += 1
----------------------------------------------------------------------


#Gera um programa gráfico. Os dois blocos de codigos abaixo devem estar em files.py separados
from Tkinter import*
class Hello(Frame):
     def __init__(self, master=None):
         Frame.__init__(self, master)
         self.pack( )
         self.make_widgets( )
         
     def make_widgets(self):
         widget = Button(self, text='Hello world', command=self.quit)
         widget.pack(side=BOTTOM)

if __name__ == '__main__': Hello( ).mainloop( )


from Dialog  import*
from Tkinter import*
from hello   import*
class HelloGoodbye(Hello):
    def really_quit(self):
        Hello.quit(self)
        
    def quit(self):
        ans = Dialog(self,
                 title = 'Verify exit',
                 text  = "Encerrar o programa.",
                 bitmap= 'question',
                 strings= ('Yes','No'), default = 1)
        if ans.num == 0:
           self.really_quit( )
           
    def make_widgets(self):
        Hello.make_widgets(self)
        extra = Button(self, text ='Eu computo', command=self.really_quit)
        extra.pack(side=TOP)
        
if __name__ == '__main__': HelloGoodbye( ).mainloop( )
-----------------------------------------------------------------------------

