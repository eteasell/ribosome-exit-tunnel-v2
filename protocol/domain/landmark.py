
class Landmark:
    def __init__(self, position: int, residue: str, name: str = None):
        self.name = name
        self.position = position
        self.residue = residue
        
    def __str__(self):
        return f"{self.residue}{self.position}"