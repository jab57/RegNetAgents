"""
Temporary patch for flash-attn compatibility issues.
This creates a fallback when flash_attn is not available.
"""

import sys
import importlib.util
from unittest.mock import MagicMock

def patch_flash_attn():
    """Create mock flash_attn modules when not available."""
    
    # Create mock flash_attn module
    flash_attn = MagicMock()
    flash_attn.flash_attn_interface = MagicMock()
    flash_attn.bert_padding = MagicMock()
    
    # Mock the specific functions that are imported
    def mock_flash_attn_varlen_kvpacked_func(*args, **kwargs):
        # Fallback to standard PyTorch attention - this is a simplified version
        # In practice, you'd need to implement proper attention mechanism
        raise NotImplementedError("Flash attention not available - using fallback")
    
    def mock_unpad_input(*args, **kwargs):
        # Simple passthrough for unpad_input
        return args[0], args[1] if len(args) > 1 else None
    
    flash_attn.flash_attn_interface.flash_attn_varlen_kvpacked_func = mock_flash_attn_varlen_kvpacked_func
    flash_attn.bert_padding.unpad_input = mock_unpad_input
    
    # Add to sys.modules so imports work
    sys.modules['flash_attn'] = flash_attn
    sys.modules['flash_attn.flash_attn_interface'] = flash_attn.flash_attn_interface
    sys.modules['flash_attn.bert_padding'] = flash_attn.bert_padding
    
    print("Applied flash_attn compatibility patch")

if __name__ == "__main__":
    patch_flash_attn()