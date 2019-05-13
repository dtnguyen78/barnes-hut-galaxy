import sys


def extract_data(file_name):
    with open(file_name, 'r') as f:
        content = f.readlines()
        content = [x.strip() for x in content[1:]]

    to_return = []
    for line in content:
        chunks = line.split(',')

        """
            os
            <<outer_idx
            <<","<<i                <-- [1]
            <<","<<t                <-- [2]
            <<","<<bodies[i].x()    <-- [3]
            <<","<<bodies[i].y()    <-- [4]
            <<","<<bodies[i].vx()   <-- [5]
            <<","<<bodies[i].vy()   <-- [6]
            <<","<<bodies[i].m()
            <<endl;
        """
        to_return.append((chunks[1], chunks[2], chunks[3], chunks[4], chunks[5], chunks[6]))
    return to_return


def compare_print(serial, barnes_hut):
    serial_data = extract_data(serial)
    bh_data = extract_data(barnes_hut)

    # both data arrays contain: body_id, time-step, x, y, vx, vy
    # merge these two and first print some stats
    for i in xrange(len(serial_data)):
        l = serial_data[i]
        r = bh_data[i]
        print "body: {} | id: {} | x_diff: {} | y_diff: {} | vx_diff: {} | vy_diff: {}" \
            .format(l[0], l[1], float(l[2]) - float(r[2]), float(l[3]) - float(r[3]),
                    float(l[4]) - float(r[4]), float(l[5]) - float(r[5]))

        # print "body: {}   id: {}   x_serial: {}  x_bh: {}  y_serial: {}  y_bh: {}  vx_serial: {}  vx_bh: {}  " \
        #       "vy_serial: {}  vy_bh: {}".format(l[0], l[1], l[2], r[2], l[3], r[3], l[4], r[4], l[5], r[5])


if __name__ == "__main__":
    if len(sys.argv) == 3:
        serial = sys.argv[1]
        barnes_hut = sys.argv[2]
        compare_print(serial, barnes_hut)
