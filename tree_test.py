from hmm.inheritance import Tree

test = Tree()

test.broaden(None, 3)

test.print_tree()

for i in range(3):
#
	test.broaden([i], 3)
#

test.print_tree()

for i in range(3):
#
	test.add_data([i], i)
	
	for j in range(3):
	#
		test.add_data([i,j], i+j)
	#
#

test.print_tree()

test.deepen([1,1], 3)

test.print_tree()

test.move_node([1,1,0,0], [1,1])

test.print_tree()

test.move_node([1,1,1,0], [1,1])

test.print_tree()

test.move_data([2,2], 0, [1,1,2])

test.print_tree()

test.deepen([1,2], 3)

test.delete_node([1,1])

test.print_tree()
