define(['jquery'], function($) {

    var $body, api = {
        debug: false
    };

    api.log = function() {
        if (this.debug) {
            console.log('[API] ' + arguments[0], Array.prototype.splice.call(arguments, 1));
        }
    };

    var Result = function (id, data) {
        this.id = id;
        this.states = data.hasOwnProperty('state') ? data['state'] : [];
        this.counts = data.hasOwnProperty('counter') ? data['counter'] : [];

        this.hasState = function(key) {
            return this.states.hasOwnProperty(key);
        };
        this.getState = function(key) {
            if (this.hasState(key)) {
                return this.states[key];
            }
            throw 'Undefined state "' + key + '"';
        };
        this.hasCount = function(key) {
            return this.counts.hasOwnProperty(key);
        };
        this.getCount = function(key) {
            if (this.hasCount(key)) {
                return this.counts[key];
            }
            throw 'Undefined count "' + key + '"';
        };
    };
    Result.prototype.getId = function() {
        return this.id;
    };
    Result.prototype.isPinned = function() {
        return this.getState('is_pinned');
    };
    Result.prototype.isBanned = function() {
        return this.getState('is_banned');
    };
    Result.prototype.isTrashed = function() {
        return this.getState('is_trashed');
    };
    Result.prototype.isAdded = function() {
        return this.getState('is_added');
    };
    Result.prototype.getPinCount = function() {
        return this.getCount('pin');
    };
    Result.prototype.getBanCount = function() {
        return this.getCount('ban');
    };
    Result.prototype.getTrashCount = function() {
        return this.getCount('trash');
    };
    Result.prototype.getLibraryCount = function() {
        return this.getCount('library');
    };

    /**
     * Toggles the pinned state of the paper.
     *
     * @param id        (int) The paper id/pk
     * @param source    (string) The source page
     */
    api.pin = function(id, source) {
        source = source || '';
        $.ajax({
            type: 'POST',
            url: '/user/paper/pin',
            data: {'id': id, 'source': source}
        })
        .done(function(data) {
            api.log('pin success', data);
            $body.trigger('etalia.publication.pin', new Result(id, data));
        })
        .fail(function(xrh, status, error) {
            api.log('pin failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Toggles the banned state of the paper.
     *
     * @param id        (int) The paper id/pk
     * @param source    (string) The source page
     */
    api.ban = function(id, source) {
        source = source || '';
        $.ajax({
            type: 'POST',
            url: '/user/paper/ban',
            data: {'id': id, 'source': source}
        })
        .done(function(data) {
            api.log('ban success', data);
            $body.trigger('etalia.publication.ban', new Result(id, data));
        })
        .fail(function(xrh, status, error) {
            api.log('ban failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Adds the paper to the user library.
     *
     * @param id (int) The paper id/pk
     * @param failureCallback (closure)
     */
    api.add = function(id, failureCallback) {
        $.ajax({
            type: 'POST',
            url: '/user/paper/add',
            data: {'id': id}
        })
        .done(function(data) {
            if (data.hasOwnProperty('success') && !data['success']) {
                if (failureCallback) {
                    failureCallback();
                }
                api.log('add failure', data);
                return;
            }
            api.log('add success', data);
            $body.trigger('etalia.publication.add', new Result(id, data));
        })
        .fail(function(xrh, status, error) {
            if (failureCallback) {
                failureCallback();
            }
            api.log('add failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Trashes the paper from the user library.
     *
     * @param id (int) The paper id/pk
     * @param failureCallback (closure)
     */
    api.trash = function(id, failureCallback) {
        $.ajax({
            type: 'POST',
            url: '/user/paper/trash',
            data: {'id': id}
        })
        .done(function(data) {
            if (data.hasOwnProperty('success') && !data['success']) {
                if (failureCallback) {
                    failureCallback();
                }
                api.log('trash failure', data);
                return;
            }
            api.log('trash success', data);
            $body.trigger('etalia.publication.trash', new Result(id, data));
        })
        .fail(function(xrh, status, error) {
            if (failureCallback) {
                failureCallback();
            }
            api.log('trash failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Restores the paper into the user library.
     *
     * @param id (int) The paper id/pk
     * @param failureCallback (closure)
     */
    api.restore = function(id, failureCallback) {
        $.ajax({
            type: 'POST',
            url: '/user/paper/restore',
            data: {'id': id}
        })
        .done(function(data) {
            if (data.hasOwnProperty('success') && !data['success']) {
                if (failureCallback) {
                    failureCallback();
                }
                api.log('restore failure', data);
                return;
            }
            api.log('restore success', data);
            $body.trigger('etalia.publication.restore', new Result(id, data));
        })
        .fail(function(xrh, status, error) {
            if (failureCallback) {
                failureCallback();
            }
            api.log('restore failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Clears the user library trash.
     */
    api.clearTrash = function() {
        $.ajax({
            type: 'POST',
            url: '/user/library/trash/empty'
        })
        .done(function(data) {
            api.log('clear trash success', data);
            $body.trigger('etalia.publication.trash-clear', new Result(undefined, data));
        })
        .fail(function(xrh, status, error) {
            api.log('clear trash failure', xrh, status, error);
        });

        return this;
    };

    $(function() {
        $body = $('body');
    });

    return api;
});
