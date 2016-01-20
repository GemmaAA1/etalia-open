define(['jquery'], function() {

    return {
        toggleClass: function($element, cssClass) {
            if ($element.hasClass(cssClass)) {
                $element.removeClass(cssClass);
                return false;
            }
            $element.addClass(cssClass);
            return true;
        }
    }

});
