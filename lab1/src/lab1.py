import math

def picar1(x):
    return x**3 / 3

def picar2(x):
    return picar1(x) + x ** 7 / 63

def picar3(x):
    return picar2(x) + (x ** 11) * (2 / 2079) + (x ** 15) / 59535
    
def picar4(x):
    return picar3(x) + (x ** 15)*(2 / 93555) + (x ** 19)*(2 / 3393495) + (x ** 19)*(2 / 2488563) + \
    (x ** 23)*(2 / 86266215) + (x ** 23)*(1 / 99411543) + (x ** 27)*(2 / 3341878155) + (x ** 31)*(1 / 109876902975)

# n - количество итераций, h - шаг, (x, y) - начальная точка
def Euler(n, h, x, y):
    answer = []
    
    for i in range(n):
        try:
            y += h * function(x, y)
            answer.append(y)   
            x += h
        except OverflowError:
            answer.append("Over")
            
    return answer # решение

def rungeKutta(n, h, x, y, alpha = 0.5):
    answer = []

    for i in range(n):
        try:
            step1 = (1 - alpha) * function(x, y)
            step2 = alpha * function(x + h / (2 * alpha), y + (h / (2 * alpha)) * function(x, y))
            step3 = step1 + step2
            y += h * step3
            answer.append(y)   
            x += h
        except OverflowError:
            answer.append("Over")
            
    return answer # решение
    
def function(x, y):
    return x ** 2 + y ** 2 # функция первой производной

def checkFormat(item):
    if type(item) == float:
        if item > 1000000:
            return '{:.2e}'.format(item)
        return '{:.2f}'.format(item)
    
    elif type(item) == int:
        return str(item)
    else:
        return item

def main():
    x_start = 0
    x_end = 3
    
    h = 1e-6 

    n = math.ceil(abs(x_end - x_start) / h) + 1 # число итераций
    
    output_step = int(n / 150) # выводим только 150 значений в таблице

    answer_euler = Euler(n, h, 0, 0)
    answer_rungeKutta = rungeKutta(n, h, 0, 0)
    
    print("------------------------------------------------------------------------------------------------------------")
    print("|         |__________________________Метод Пикара_________________________|    Метод     |      Метод      |")
    print("|    x    |               |               |               |               |    Эйлера    |   Рунге-Кутты   |")
    print("|         |   1-е прибл.  |   2-е прибл.  |   3-е прибл.  |   4-е прибл.  |              |                 |")
    print("------------------------------------------------------------------------------------------------------------")
    
    for i in range(0, n, output_step):
         print("|{:^9.2f}|{:^15.2f}|{:^15.2f}|{:^15.2f}|{:^15.2f}|{:^14s}|{:^17s}|".format(x_start, picar1(x_start), picar2(x_start), \
                                                                                     picar3(x_start), picar4(x_start), \
                                                                                     checkFormat(answer_euler[i]), checkFormat(answer_rungeKutta[i])))
         x_start += h * output_step
        
        
main()
