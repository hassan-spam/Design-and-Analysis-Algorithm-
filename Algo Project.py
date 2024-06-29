import numpy as np
import time
import matplotlib.pyplot as plt
from memory_profiler import profile
from functools import reduce
from scipy.optimize import fsolve

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def calculate_vector(x1, y1, x2, y2):
    return (x2 - x1, y2 - y1)

def calculate_cross_product(v1, v2):
    return v1[0] * v2[1] - v1[1] * v2[0]

def check_intersection_cross(segment1, segment2):
    x1, y1, x2, y2 = segment1
    x3, y3, x4, y4 = segment2
    
    vector1 = calculate_vector(x1, y1, x2, y2)
    vector2 = calculate_vector(x3, y3, x4, y4)
    
    cross_product = calculate_cross_product(vector1, vector2)
    
    if cross_product == 0:
        return True
    
    vector3 = calculate_vector(x1, y1, x3, y3)
    vector4 = calculate_vector(x1, y1, x4, y4)
    
    if (
        (cross_product > 0 and calculate_cross_product(vector3, vector4) < 0) or
        (cross_product < 0 and calculate_cross_product(vector3, vector4) > 0)
    ):
        return True

    return False

def calculate_slope(x1, y1, x2, y2):
    if x1 == x2:
        return None
    return (y2 - y1) / (x2 - x1)

def check_intersection_slope(segment1, segment2):
    x1, y1, x2, y2 = segment1
    x3, y3, x4, y4 = segment2
    
    slope1 = calculate_slope(x1, y1, x2, y2)
    slope2 = calculate_slope(x3, y3, x4, y4)
    
    if slope1 is None and slope2 is None:
        if x1 == x3 and (min(y1, y2) <= max(y3, y4)) and (min(y3, y4) <= max(y1, y2)):
            return True
        else:
            return False
    elif slope1 is None or slope2 is None:
        return False
    
    if slope1 == slope2:
        if (x1 == x3) and (y1 == y3) and (x2 == x4) and (y2 == y4):
            return True
        if (min(x1, x2) <= max(x3, x4)) and (min(x3, x4) <= max(x1, x2)) and \
           (min(y1, y2) <= max(y3, y4)) and (min(y3, y4) <= max(y1, y2)):
            return True

    return False

def plot_segments(segment1, segment2, intersection):
    x1, y1, x2, y2 = segment1
    x3, y3, x4, y4 = segment2

    plt.plot([x1, x2], [y1, y2], 'b-', label='Segment 1')
    plt.plot([x3, x4], [y3, y4], 'g-', label='Segment 2')

    if intersection:
        intersection_x, intersection_y = intersection
        plt.scatter(intersection_x, intersection_y, color='red', marker='o', label='Intersection')
    
    plt.legend()
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Line Segment Intersection')
    plt.grid(True)
    plt.show()

def find_intersection_point(segment1, segment2):
    x1, y1, x2, y2 = segment1
    x3, y3, x4, y4 = segment2

    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

    if denominator == 0:
        return None

    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denominator
    intersection_x = x1 + t * (x2 - x1)
    intersection_y = y1 + t * (y2 - y1)

    if (
        min(x1, x2) <= intersection_x <= max(x1, x2) and
        min(y1, y2) <= intersection_y <= max(y1, y2) and
        min(x3, x4) <= intersection_x <= max(x3, x4) and
        min(y3, y4) <= intersection_y <= max(y3, y4)
    ):
        return intersection_x, intersection_y
    else:
        return None

def find_intersection_numeric(segment1, segment2):
    def equations(vars):
        x, y = vars
        eq1 = (y - segment1[1]) - ((segment1[3] - segment1[1]) / (segment1[2] - segment1[0])) * (x - segment1[0])
        eq2 = (y - segment2[1]) - ((segment2[3] - segment2[1]) / (segment2[2] - segment2[0])) * (x - segment2[0])
        return [eq1, eq2]

    initial_guess = [(segment1[0] + segment2[0]) / 2, (segment1[1] + segment2[1]) / 2]

    intersection_point = fsolve(equations, initial_guess)

    return tuple(intersection_point)


def orientation(p, q, r):
    val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
    if val == 0:
        return 0  # collinear
    return 1 if val > 0 else 2  # clock or counterclockwise

def monotone_chain_convex_hull(points):
    n = len(points)
    
    points.sort(key=lambda point: (point.x, point.y))


    lower_hull = []
    for p in points:
        while len(lower_hull) >= 2 and orientation(lower_hull[-2], lower_hull[-1], p) != 2:
            lower_hull.pop()
        lower_hull.append(p)


    upper_hull = []
    for p in reversed(points):
        while len(upper_hull) >= 2 and orientation(upper_hull[-2], upper_hull[-1], p) != 2:
            upper_hull.pop()
        upper_hull.append(p)

   
    convex_hull = lower_hull[:-1] + upper_hull[:-1]
    
def on_segment(p, q, r):
    return (q.x <= max(p.x, r.x) and q.x >= min(p.x, r.x) and
            q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y))

def do_intersect(p1, q1, p2, q2):
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    if o1 != o2 and o3 != o4:
        return True

    if o1 == 0 and on_segment(p1, p2, q1):
        return True
    if o2 == 0 and on_segment(p1, q2, q1):
        return True
    if o3 == 0 and on_segment(p2, p1, q2):
        return True
    if o4 == 0 and on_segment(p2, q1, q2):
        return True

    return False

def brute_force_convex_hull(points):
    n = len(points)
    convex_hull = []

    for i in range(n):
        for j in range(i+1, n):
            valid = True
            for k in range(n):
                if k != i and k != j:
                    if orientation(points[i], points[j], points[k]) == 0 and on_segment(points[i], points[k], points[j]):
                        valid = False
                        break

            if valid:
                if not any(do_intersect(points[i], points[j], convex_hull[p], convex_hull[(p+1) % len(convex_hull)]) for p in range(len(convex_hull))):
                    convex_hull.append(points[i])
                    convex_hull.append(points[j])

    return convex_hull
    
def direction(p, q, r):
    return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)

def distance_sq(p, q):
    return (p.x - q.x) ** 2 + (p.y - q.y) ** 2

def jarvis_march(points):
    
    a = min(points, key=lambda point: point.x)
    index = points.index(a)

    
    l = index
    result = []
    result.append(a)
    while True:
        q = (l + 1) % len(points)
        for i in range(len(points)):
            if i == l:
                continue
            
            d = direction(points[l], points[i], points[q])
            if d > 0 or (d == 0 and distance_sq(points[i], points[l]) > distance_sq(points[q], points[l])):
                q = i
        l = q
        if l == index:
            break
        result.append(points[q])

    return result

def convex_hull_graham(points):
    TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

    def cmp(a, b):
        return (a > b) - (a < b)

    def turn(p, q, r):
        return cmp((q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]), 0)

    def _keep_left(hull, r):
        while len(hull) > 1 and turn(hull[-2], hull[-1], r) != TURN_LEFT:
            hull.pop()
        if not len(hull) or hull[-1] != r:
            hull.append(r)
        return hull

    points = sorted(points)
    l = reduce(_keep_left, points, [])
    u = reduce(_keep_left, reversed(points), [])
    convex_hull = l + [u[i] for i in range(1, len(u) - 1)]
    display_convex_hull(points, convex_hull)
    return convex_hull

def quickhull(points):
    def get_distance(p1, p2, p):
        return (p[1] - p1[1]) * (p2[0] - p1[0]) - (p2[1] - p1[1]) * (p[0] - p1[0])

    def find_hull(p1, p2, points, side, hull):
        max_dist = -1
        farthest_point = None
        
        for p in points:
            dist = get_distance(p1, p2, p)
            if dist * side > 0:
                distance = abs(dist)
                if distance > max_dist:
                    max_dist = distance
                    farthest_point = p
        
        if farthest_point:
            index = points.index(farthest_point)
            hull.append(farthest_point)
            points.remove(farthest_point)
            find_hull(p1, farthest_point, points, -get_distance(p1, farthest_point, p2), hull)
            find_hull(farthest_point, p2, points, -get_distance(farthest_point, p2, p1), hull)

    if len(points) < 3:
        return "Need at least 3 points for a convex hull"
    
    leftmost = min(points, key=lambda x: x[0])
    rightmost = max(points, key=lambda x: x[0])

    convex_hull = [leftmost, rightmost]
    points.remove(leftmost)
    points.remove(rightmost)

    left_side = []
    right_side = []
    
    for point in points:
        side = get_distance(leftmost, rightmost, point)
        if side > 0:
            left_side.append(point)
        elif side < 0:
            right_side.append(point)

    find_hull(leftmost, rightmost, left_side, -1, convex_hull)
    find_hull(rightmost, leftmost, right_side, -1, convex_hull)

    
    return convex_hull


def display_convex_hull(points, convex_hull):
    # Plotting the points and convex hull
    x, y = zip(*[(point.x, point.y) for point in points])
    hull_x, hull_y = zip(*[(point.x, point.y) for point in convex_hull])

    plt.figure(figsize=(8, 8))
    plt.scatter(x, y, color='blue', label='Points')
    plt.plot(hull_x + (hull_x[0],), hull_y + (hull_y[0],), color='red', label='Convex Hull')
    plt.title('Convex Hull')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.legend()
    plt.grid(True)
    plt.show()
    return convex_hull

def get_segment_input():
    segment1 = tuple(map(int, input("Enter segment 1 as x1, y1, x2, y2: ").split(',')))
    segment2 = tuple(map(int, input("Enter segment 2 as x1, y1, x2, y2: ").split(',')))
    return segment1, segment2

@profile
def measure_time_and_space_numeric():
    start_time = time.process_time()
    segment1, segment2 = get_segment_input()
    intersection_numeric = find_intersection_numeric(segment1, segment2)
    print("\nUsing the numerical method:")
    if intersection_numeric:
        print("Line segments do not intersect.")
    else:
        print("Line segments  intersect.")
        
    plot_segments(segment1, segment2, intersection_numeric)
    end_time = time.process_time()

    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time:.6f} seconds")
    
@profile    
def measure_time_and_space_slope():
    start_time = time.process_time()
    segment1, segment2 = get_segment_input()

    intersection_slope = check_intersection_slope(segment1, segment2)
    print("Using the slope method:")
    if intersection_slope:
        print("Line segments intersect.")
    else:
        print("Line segments do not intersect.")

    intersection_point = find_intersection_point(segment1, segment2)
    plot_segments(segment1, segment2, intersection_point)
    end_time = time.process_time()

    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time:.6f} seconds")


@profile
def measure_time_and_space_cross():
    start_time = time.process_time()
    segment1, segment2 = get_segment_input()

    # Your existing code
    check_intersection_cross(segment1, segment2)
    intersection_cross = check_intersection_cross(segment1, segment2)
    print("\nUsing the cross product method:")
    if intersection_cross:
        print("Line segments do not intersect.")
    else:
        print("Line segments intersect.")

    intersection_point = find_intersection_point(segment1, segment2)
    plot_segments(segment1, segment2, intersection_point)

    end_time = time.process_time()

    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time:.6f} seconds")


def get_points_from_user():
    num_points = int(input("Enter the number of points: "))
    points = []

    for i in range(num_points):
        x, y = map(int, input(f"Enter x and y coordinates for point {i + 1}: ").split())
        point = Point(x, y)
        points.append(point)

    return points

def measure_time_and_space_monotone():
    start_time = time.process_time()
    points = get_points_from_user()
    convex_hull = monotone_chain_convex_hull(points)
    display_convex_hull(points, convex_hull)
    end_time = time.process_time()

    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time:.6f} seconds")


def measure_time_and_space_jarvis():
    start_time = time.process_time()
    points = get_points_from_user()
    convex_hull = jarvis_march(points)
    display_convex_hull(points, convex_hull)
    end_time = time.process_time()

    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time:.6f} seconds")

def measure_time_and_space_graham():
    start_time = time.process_time()
    points = get_points_from_user()
    convex_hull = convex_hull_graham(points)
    display_convex_hull(points, convex_hull)
    end_time = time.process_time()

    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time:.6f} seconds")


def measure_time_and_space_brute():
    start_time = time.process_time()
    points = get_points_from_user()
    convex_hull = brute_force_convex_hull(points)
    display_convex_hull(points, convex_hull)
    end_time = time.process_time()

    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time:.6f} seconds")
    

def measure_time_and_space_quick():
    start_time = time.process_time()
    points = get_points_from_user()
    convex_hull = quickhull(points)
    display_convex_hull(points, convex_hull)
    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time:.6f} seconds")
    


def convex_hull():
    while True:
        print("Please select the function you want to measure time and space for:")
        print("1. brute Force Method")
        print("2. quick Hull Method")
        print("3. monotone Method")
        print("4. graham Method")
        print("5. jarvis Method")
        print("0. Exit")
        user_choice = int(input("Enter your choice (1-5): "))

        if user_choice == 1:
            measure_time_and_space_brute()
        elif user_choice == 2:
            measure_time_and_space_quick()
        elif user_choice == 3:
            measure_time_and_space_monotone()
        elif user_choice == 4:
            measure_time_and_space_graham()
        elif user_choice == 5:
            measure_time_and_space_jarvis()
        else:
            print("Invalid choice. Please select a valid option (1-5).")
        
def intersection():
    while True:
        print("\nChoose a function to call:")
        print("1.  Cross Method")
        print("2.  Slope Method")
        print("3.  Numeric Method")
        print("0. Exit")
        choice = input("Enter your choice (0 to exit, 1-3 to call functions): ")

        if choice == '0':
            print("Exiting the program.")
            break
        elif choice == '1':
            measure_time_and_space_cross()
        elif choice == '2':
            measure_time_and_space_slope()
        elif choice == '3':
            measure_time_and_space_numeric()
        else:
            print("Invalid choice. Please enter 0, 1, 2, or 3.")

while True:
        print("\nChoose an option:")
        print("1. Intersection")
        print("2. Convex Hull")
        print("0. Exit")
        option = int(input("Enter your choice (0 to exit, 1 or 2): "))

        if option == 0:
            print("Exiting the program.")
            break
        elif option == 1:
            intersection()
        elif option == 2:
            convex_hull()
        else:
            print("Invalid choice. Please enter 0, 1, or 2.")
